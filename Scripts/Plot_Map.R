isolates <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(reference_strain ==1 ) 
  dplyr::filter(isotype %in% names(strain_islands) & !(isotype %in% hm_hi1$isotype)) %>%
  dplyr::mutate(altitude = as.numeric(NA)) %>%
  dplyr::mutate(longitude = as.numeric(longitude),
                substrate_temperature = substrate_temp,
                ambient_temperature = ambient_temp)