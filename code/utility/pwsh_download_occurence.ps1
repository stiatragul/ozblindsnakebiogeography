# In the data folder run this script

Invoke-WebRequest -Uri "https://doi.ala.org.au/doi/d92678b1-ad2d-437b-9457-9f52737ba003/download" -Outfile "2021_ALA_blindsnake_occurence_data.zip"

Expand-Archive .\2021_ALA_blindsnake_occurence_data.zip

Remove-Item .\2021_ALA_blindsnake_occurence_data.zip

# Copy-Item "./2021_ALA_blindsnake_occurence_data/records*.csv" -Destination "./records.csv"