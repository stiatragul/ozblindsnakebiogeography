@ECHO OFF
ECHO Two way mirror between local_repo and OneDrive for large files *.ai *.psd *.pdf *.png *.jpg *.png *.xlsx *.xls *.docx *.qgz 

robocopy "C:\Users\ST\Documents\repo\blindsnakebiogeography\\" "C:\Users\ST\OneDrive - Australian National University\repo_onedrive\blindsnakebiogeography_onedrive\\" *.ai *.psd *.pdf *.png *.jpg *.png *.xlsx *.xls *.docx *.qgz *wwf_terr* *map* terr_biome* *.RData *.dbf *.prj *.shp *.shx *.qpj *.ait *.zip *.jld2 *.svg /S /XO /R:1 /W:1 /NDL /XJD /v /LOG+:"C:\Users\ST\Documents\repo\blindsnakebiogeography:src.log"

robocopy "C:\Users\ST\OneDrive - Australian National University\repo_onedrive\blindsnakebiogeography_onedrive\\" "C:\Users\ST\OneDrive - Australian National University\repo_onedrive\blindsnakebiogeography_onedrive\\" *.ai *.psd *.pdf *.png *.jpg *.png *.xlsx *.xls *.docx *.qgz *wwf_terr* *map* terr_biome* *.RData *.dbf *.prj *.shp *.shx *.qpj *.ait *.zip *.jld2 *.svg /S /XO /R:1 /W:1 /NDL /XJD /v /LOG+:"C:\Users\ST\OneDrive - Australian National University\repo_onedrive\blindsnakebiogeography_onedrive\\src.log"

echo Sync Complete