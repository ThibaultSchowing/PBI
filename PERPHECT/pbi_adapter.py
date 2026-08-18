"""
Adapter bridging PBI-Scope data to PERPHECT model format.

Converts PBI-Scope's string-based phage/host IDs and raw sequences into
the integer-indexed, zero-padded, one-hot-encoded format that PERPHECT expects.

Usage:
    from pbi import quick_connect
    from pbi_adapter import PBIAdapter

    retriever = quick_connect()
    adapter = PBIAdapter(retriever)

    # Get DataFrames in PERPHECT format
    couples_df, bacteria_df, phages_df = adapter.to_perphect_dataframes(
        positive_pairs, negative_pairs
    )

    # Or use the TF generator for memory-efficient training
    gen = adapter.create_tf_generator(couples, labels, batch_size=16)
"""

import math
import logging
from typing import Optional, Tuple, Dict, List, Generator

import numpy as np
import pandas as pd

from transforms import translate_sequence_onehot

logger = logging.getLogger(__name__)

# PERPHECT thresholds (bp)
BACTERIUM_THRESHOLD = 7_000_000
PHAGE_THRESHOLD = 200_000
BACTERIUM_MIN_LENGTH = 150_000
PHAGE_MIN_LENGTH = 1_500

# Interaction types that map to label=0 (negative)
NEGATIVE_INTERACTIONS = {"no interaction", "none", "negative", "non-interacting"}


class PBIAdapter:
    """
    Bridges PBI-Scope data to PERPHECT model format.

    Handles ID mapping, sequence padding, one-hot encoding, and provides
    both DataFrame-based and generator-based data access.

    Args:
        retriever: PBI-Scope SequenceRetriever instance.
        bacterium_threshold: Max bacteria sequence length after padding (default 7M).
        phage_threshold: Max phage sequence length after padding (default 200K).
        bacterium_min_length: Minimum bacteria length to keep (default 150K).
        phage_min_length: Minimum phage length to keep (default 1.5K).
    """

    def __init__(
        self,
        retriever,
        bacterium_threshold: int = BACTERIUM_THRESHOLD,
        phage_threshold: int = PHAGE_THRESHOLD,
        bacterium_min_length: int = BACTERIUM_MIN_LENGTH,
        phage_min_length: int = PHAGE_MIN_LENGTH,
    ):
        self.retriever = retriever
        self.bacterium_threshold = bacterium_threshold
        self.phage_threshold = phage_threshold
        self.bacterium_min_length = bacterium_min_length
        self.phage_min_length = phage_min_length

        # ID mapping: PBI string ID -> PERPHECT integer ID
        self._host_id_map: Dict[str, int] = {}
        self._phage_id_map: Dict[str, int] = {}
        self._next_host_id = 0
        self._next_phage_id = 0

        # Sequence caches
        self._host_sequences: Dict[str, str] = {}
        self._phage_sequences: Dict[str, str] = {}

    def _map_host_id(self, pbi_host_id: str) -> int:
        """Map a PBI-Scope Host_ID string to a PERPHECT integer ID."""
        if pbi_host_id not in self._host_id_map:
            self._host_id_map[pbi_host_id] = self._next_host_id
            self._next_host_id += 1
        return self._host_id_map[pbi_host_id]

    def _map_phage_id(self, pbi_phage_id: str) -> int:
        """Map a PBI-Scope Phage_ID string to a PERPHECT integer ID."""
        if pbi_phage_id not in self._phage_id_map:
            self._phage_id_map[pbi_phage_id] = self._next_phage_id
            self._next_phage_id += 1
        return self._phage_id_map[pbi_phage_id]

    def _fetch_host_sequence(self, host_id: str) -> Optional[str]:
        """Fetch and cache a host sequence from PBI-Scope."""
        if host_id in self._host_sequences:
            return self._host_sequences[host_id]

        try:
            seq = self.retriever.get_host_sequence(host_id, contig_mode="concat")
            if seq and len(seq) >= self.bacterium_min_length:
                self._host_sequences[host_id] = seq
                return seq
            else:
                logger.debug(
                    f"Host {host_id} too short ({len(seq) if seq else 0} bp), "
                    f"minimum is {self.bacterium_min_length}"
                )
                return None
        except Exception as e:
            logger.warning(f"Failed to fetch host sequence for {host_id}: {e}")
            return None

    def _fetch_phage_sequence(self, phage_id: str) -> Optional[str]:
        """Fetch and cache a phage sequence from PBI-Scope."""
        if phage_id in self._phage_sequences:
            return self._phage_sequences[phage_id]

        try:
            seq = self.retriever.get_phage_sequence(phage_id)
            if seq and len(seq) >= self.phage_min_length:
                self._phage_sequences[phage_id] = seq
                return seq
            else:
                logger.debug(
                    f"Phage {phage_id} too short ({len(seq) if seq else 0} bp), "
                    f"minimum is {self.phage_min_length}"
                )
                return None
        except Exception as e:
            logger.warning(f"Failed to fetch phage sequence for {phage_id}: {e}")
            return None

    def _pad_and_encode(self, sequence: str, max_length: int) -> np.ndarray:
        """
        Zero-pad (or truncate) a sequence and convert to one-hot encoding.

        Args:
            sequence: Raw DNA sequence string.
            max_length: Target length after padding.

        Returns:
            numpy array of shape (max_length, 4), dtype uint8.
        """
        # Truncate if too long
        truncated = sequence[:max_length]
        # One-hot encode
        onehot = translate_sequence_onehot(truncated)
        # Zero-pad if needed
        if onehot.shape[0] < max_length:
            padding = np.zeros((max_length - onehot.shape[0], 4), dtype=np.uint8)
            onehot = np.concatenate([onehot, padding], axis=0)
        return onehot

    # ------------------------------------------------------------------
    # DataFrame conversion
    # ------------------------------------------------------------------

    def to_perphect_dataframes(
        self,
        positive_pairs: pd.DataFrame,
        negative_pairs: Optional[pd.DataFrame] = None,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Convert PBI-Scope pairs to PERPHECT-format DataFrames.

        Args:
            positive_pairs: DataFrame with Phage_ID, Host_ID, and optionally
                           Phage_Sequence, Host_Sequence columns.
            negative_pairs: Optional DataFrame with same columns plus Label=0.

        Returns:
            Tuple of (couples_df, bacteria_df, phages_df) matching PERPHECT
            CSV format:
            - couples_df: columns [id, bacterium_id, phage_id, interaction_type]
            - bacteria_df: index=bacterium_id, column=bacterium_sequence (one-hot)
            - phages_df: index=phage_id, column=phage_sequence (one-hot)
        """
        records = []
        bacteria_seqs = {}
        phage_seqs = {}

        # Process positive pairs
        for _, row in positive_pairs.iterrows():
            phage_id_str = row["Phage_ID"]
            host_id_str = row["Host_ID"]

            phage_id_int = self._map_phage_id(phage_id_str)
            host_id_int = self._map_host_id(host_id_str)

            # Fetch sequences if not provided
            phage_seq = row.get("Phage_Sequence")
            host_seq = row.get("Host_Sequence")

            if phage_seq is None or (isinstance(phage_seq, float) and pd.isna(phage_seq)):
                phage_seq = self._fetch_phage_sequence(phage_id_str)
            if host_seq is None or (isinstance(host_seq, float) and pd.isna(host_seq)):
                host_seq = self._fetch_host_sequence(host_id_str)

            if phage_seq is None or host_seq is None:
                continue

            # Filter by minimum length
            if len(host_seq) < self.bacterium_min_length:
                continue
            if len(phage_seq) < self.phage_min_length:
                continue

            # Store sequences (will be one-hot encoded later)
            bacteria_seqs[host_id_int] = host_seq
            phage_seqs[phage_id_int] = phage_seq

            records.append({
                "bacterium_id": host_id_int,
                "phage_id": phage_id_int,
                "interaction_type": 1,
            })

        # Process negative pairs
        if negative_pairs is not None:
            for _, row in negative_pairs.iterrows():
                phage_id_str = row["Phage_ID"]
                host_id_str = row["Host_ID"]

                phage_id_int = self._map_phage_id(phage_id_str)
                host_id_int = self._map_host_id(host_id_str)

                phage_seq = row.get("Phage_Sequence")
                host_seq = row.get("Host_Sequence")

                if phage_seq is None or (isinstance(phage_seq, float) and pd.isna(phage_seq)):
                    phage_seq = self._fetch_phage_sequence(phage_id_str)
                if host_seq is None or (isinstance(host_seq, float) and pd.isna(host_seq)):
                    host_seq = self._fetch_host_sequence(host_id_str)

                if phage_seq is None or host_seq is None:
                    continue

                if len(host_seq) < self.bacterium_min_length:
                    continue
                if len(phage_seq) < self.phage_min_length:
                    continue

                bacteria_seqs[host_id_int] = host_seq
                phage_seqs[phage_id_int] = phage_seq

                records.append({
                    "bacterium_id": host_id_int,
                    "phage_id": phage_id_int,
                    "interaction_type": 0,
                })

        # Build couples DataFrame
        couples_df = pd.DataFrame(records)
        if len(couples_df) > 0:
            couples_df.insert(0, "id", range(len(couples_df)))
        else:
            couples_df = pd.DataFrame(columns=["id", "bacterium_id", "phage_id", "interaction_type"])

        # Build bacteria DataFrame with one-hot encoded sequences
        if bacteria_seqs:
            bacteria_data = []
            for bacterium_id, seq in bacteria_seqs.items():
                bacteria_data.append({
                    "bacterium_id": bacterium_id,
                    "bacterium_sequence": self._pad_and_encode(seq, self.bacterium_threshold),
                })
            bacteria_df = pd.DataFrame(bacteria_data).set_index("bacterium_id")
        else:
            bacteria_df = pd.DataFrame(columns=["bacterium_id", "bacterium_sequence"]).set_index("bacterium_id")

        # Build phages DataFrame with one-hot encoded sequences
        if phage_seqs:
            phage_data = []
            for phage_id, seq in phage_seqs.items():
                phage_data.append({
                    "phage_id": phage_id,
                    "phage_sequence": self._pad_and_encode(seq, self.phage_threshold),
                })
            phages_df = pd.DataFrame(phage_data).set_index("phage_id")
        else:
            phages_df = pd.DataFrame(columns=["phage_id", "phage_sequence"]).set_index("phage_id")

        logger.info(
            f"Built PERPHECT dataframes: {len(couples_df)} couples, "
            f"{len(bacteria_df)} bacteria, {len(phages_df)} phages"
        )

        return couples_df, bacteria_df, phages_df

    # ------------------------------------------------------------------
    # TensorFlow generator (matching PERPHECT's generator signature)
    # ------------------------------------------------------------------

    def create_tf_generator(
        self,
        couples: np.ndarray,
        labels: np.ndarray,
        min_index: int = 0,
        max_index: Optional[int] = None,
        shuffle: bool = True,
        batch_size: int = 16,
    ) -> Generator:
        """
        Create a TensorFlow-compatible generator that yields batches.

        This generator matches the signature of PERPHECT's original generator()
        function, yielding ([bacterium_samples, phage_samples], targets).

        Args:
            couples: numpy array of shape (N, 2) with [bacterium_id, phage_id]
                     as integer IDs.
            labels: numpy array of shape (N,) with interaction labels (0 or 1).
            min_index: Start index for this generator.
            max_index: End index (exclusive). Defaults to len(couples).
            shuffle: Whether to randomly sample within the range.
            batch_size: Number of samples per batch.

        Yields:
            Tuple of ([bacterium_batch, phage_batch], target_batch) where
            bacterium_batch has shape (batch_size, BACTERIUM_THRESHOLD, 4),
            phage_batch has shape (batch_size, PHAGE_THRESHOLD, 4).
        """
        if max_index is None:
            max_index = len(couples)

        # Build reverse maps: integer ID -> PBI string ID
        reverse_host_map = {v: k for k, v in self._host_id_map.items()}
        reverse_phage_map = {v: k for k, v in self._phage_id_map.items()}

        while True:
            if shuffle:
                indices = np.random.randint(min_index, max_index, size=batch_size)
            else:
                indices = np.arange(min_index, min(min_index + batch_size, max_index))
                min_index += len(indices)
                if min_index >= max_index:
                    min_index = 0

            bacterium_samples = np.zeros(
                (len(indices), self.bacterium_threshold, 4), dtype=np.uint8
            )
            phage_samples = np.zeros(
                (len(indices), self.phage_threshold, 4), dtype=np.uint8
            )
            targets = np.zeros((len(indices),))

            for j, idx in enumerate(indices):
                host_id_int = couples[idx, 0]
                phage_id_int = couples[idx, 1]

                host_id_str = reverse_host_map.get(host_id_int)
                phage_id_str = reverse_phage_map.get(phage_id_int)

                if host_id_str is None or phage_id_str is None:
                    continue

                host_seq = self._host_sequences.get(host_id_str)
                phage_seq = self._phage_sequences.get(phage_id_str)

                if host_seq is None or phage_seq is None:
                    continue

                bacterium_samples[j] = self._pad_and_encode(
                    host_seq, self.bacterium_threshold
                )
                phage_samples[j] = self._pad_and_encode(
                    phage_seq, self.phage_threshold
                )
                targets[j] = labels[idx]

            yield (bacterium_samples, phage_samples), targets

    # ------------------------------------------------------------------
    # Helper: prepare training data from PBI-Scope
    # ------------------------------------------------------------------

    def prepare_training_data(
        self,
        positive_pairs: pd.DataFrame,
        negative_pairs: Optional[pd.DataFrame] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepare integer-indexed couples and labels arrays for the generator.

        This loads and caches all sequences, filters by length, and returns
        the couple/label arrays needed by create_tf_generator().

        Args:
            positive_pairs: DataFrame with Phage_ID, Host_ID columns.
            negative_pairs: Optional negative pairs DataFrame.

        Returns:
            Tuple of (couples_array, labels_array) where couples_array is
            shape (N, 2) with integer IDs and labels_array is shape (N,).
        """
        records = []

        # Process positives
        for _, row in positive_pairs.iterrows():
            phage_id_str = row["Phage_ID"]
            host_id_str = row["Host_ID"]

            phage_seq = self._fetch_phage_sequence(phage_id_str)
            host_seq = self._fetch_host_sequence(host_id_str)

            if phage_seq is None or host_seq is None:
                continue
            if len(host_seq) < self.bacterium_min_length:
                continue
            if len(phage_seq) < self.phage_min_length:
                continue

            phage_id_int = self._map_phage_id(phage_id_str)
            host_id_int = self._map_host_id(host_id_str)

            records.append((host_id_int, phage_id_int, 1))

        # Process negatives
        if negative_pairs is not None:
            for _, row in negative_pairs.iterrows():
                phage_id_str = row["Phage_ID"]
                host_id_str = row["Host_ID"]

                phage_seq = self._fetch_phage_sequence(phage_id_str)
                host_seq = self._fetch_host_sequence(host_id_str)

                if phage_seq is None or host_seq is None:
                    continue
                if len(host_seq) < self.bacterium_min_length:
                    continue
                if len(phage_seq) < self.phage_min_length:
                    continue

                phage_id_int = self._map_phage_id(phage_id_str)
                host_id_int = self._map_host_id(host_id_str)

                records.append((host_id_int, phage_id_int, 0))

        if not records:
            raise ValueError(
                "No valid pairs found after filtering. "
                "Check that sequences meet minimum length requirements "
                f"(bacteria >= {self.bacterium_min_length}, "
                f"phages >= {self.phage_min_length})."
            )

        records_arr = np.array(records)
        couples = records_arr[:, :2].astype(np.int64)
        labels = records_arr[:, 2].astype(np.float32)

        logger.info(
            f"Prepared {len(couples)} pairs "
            f"({int(labels.sum())} positive, {int(len(labels) - labels.sum())} negative)"
        )

        return couples, labels

    @property
    def host_id_map(self) -> Dict[str, int]:
        """Return the current host ID mapping."""
        return dict(self._host_id_map)

    @property
    def phage_id_map(self) -> Dict[str, int]:
        """Return the current phage ID mapping."""
        return dict(self._phage_id_map)

    def get_id_maps(self) -> Tuple[Dict[str, int], Dict[str, int]]:
        """Return both ID mappings as (host_map, phage_map)."""
        return self.host_id_map, self.phage_id_map

    # ------------------------------------------------------------------
    # Interaction type classification
    # ------------------------------------------------------------------

    def classify_pairs_by_interaction(
        self, all_pairs: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Classify pairs into positives and negatives based on interaction type.

        Queries the private_interactions table for interaction types. Pairs
        with interaction in NEGATIVE_INTERACTIONS become negatives (label=0);
        all others become positives (label=1). This is source-agnostic — any
        future source (PhageScope, etc.) that records interaction types in
        private_interactions will be handled automatically.

        Args:
            all_pairs: DataFrame from get_phage_host_pairs() with at least
                       Phage_ID and Host_ID columns.

        Returns:
            Tuple of (positive_pairs, negative_pairs) where negative_pairs
            has an additional 'negative_source' column set to 'private_data'.
        """
        try:
            interaction_df = self.retriever.conn.execute(
                "SELECT Phage_ID, Host_ID, LOWER(TRIM(interaction)) as interaction "
                "FROM private_interactions"
            ).fetchdf()
        except Exception as e:
            logger.warning(f"Could not query private_interactions: {e}")
            logger.info("Treating all pairs as positive")
            return all_pairs, pd.DataFrame(
                columns=list(all_pairs.columns) + ["negative_source"]
            )

        if interaction_df.empty:
            logger.info("No entries in private_interactions — all pairs treated as positive")
            return all_pairs, pd.DataFrame(
                columns=list(all_pairs.columns) + ["negative_source"]
            )

        # Build lookup: (Phage_ID, Host_ID) -> interaction type
        interaction_map = dict(
            zip(
                zip(interaction_df["Phage_ID"], interaction_df["Host_ID"]),
                interaction_df["interaction"],
            )
        )

        positives = []
        negatives = []
        for _, row in all_pairs.iterrows():
            key = (row["Phage_ID"], row["Host_ID"])
            interaction = interaction_map.get(key)

            if interaction in NEGATIVE_INTERACTIONS:
                neg_row = row.copy()
                neg_row["negative_source"] = "private_data"
                negatives.append(neg_row)
            else:
                positives.append(row)

        pos_df = (
            pd.DataFrame(positives)
            if positives
            else pd.DataFrame(columns=all_pairs.columns)
        )
        neg_df = (
            pd.DataFrame(negatives)
            if negatives
            else pd.DataFrame(columns=list(all_pairs.columns) + ["negative_source"])
        )

        logger.info(
            f"Classified {len(all_pairs)} pairs: "
            f"{len(pos_df)} positive, {len(neg_df)} negative "
            f"(from private_data)"
        )

        return pos_df, neg_df
