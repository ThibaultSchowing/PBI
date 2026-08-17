# PERPHECT 

Perfect is an old Phage generator that also contains a Phage-Host (Phage-Bacteria) interaction predictor. In order to test PBI-Scope in a real usage setting, we want here to adapt our data streaming method to serve the PERPHECT predictor training. 

The data will thus come directly from PBI and will have to be processed (padded) on the fly ! A good test opportunity ! 

Note: the phage, bacteria and interaction df should have these headers: "phage_id,phage_sequence", "bacterium_id,bacterium_sequence" and "id,bacterium_id,phage_id,interaction_type" (interaction type is 0 or 1 here, attention to adapt this well !
