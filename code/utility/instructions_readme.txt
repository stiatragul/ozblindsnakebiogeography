Before analysis

1. Download occurence data from Atlas of the Living Australia

```
powershell pwsh_download_occurence.ps1
```

2. Clean occurence data filter and omit some -- occurence_data_clean.R

3. Decide which bioregions/biome the species should be in -- designate_bioregions.R

4. Check sampling bias to see where more data can be collected using sampling_bias.R

5. Filter for species we have phylogenetic data for -- tree_data_match.R

6. Install biogeobears from Github -- install_biogeobears.R