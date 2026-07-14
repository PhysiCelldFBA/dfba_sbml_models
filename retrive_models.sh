echo "Dowloading C. beijerinckii (iCB925) model"
wget "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1507180034/2/MODEL1507180034_url.xml" --output-document models/C_iCB925_beijerinckii.xml
echo "Model stored in models/C_beijerinckii_iCB925.xml"
echo

echo "Dowaloading M. barkeri (iAF692) model"
wget "http://bigg.ucsd.edu/static/models/iAF692.xml" --output-document models/M_barkeri_iAF692.xml
echo "Model stored in models/M_barkeri_iAF692.xml"

echo "Dowaloading M. barkeri () model"
echo "Manually download the model from the Supporting information from:"
echo " - https://analyticalsciencejournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fbiot.201200266&file=biot_201200266_sm_suppinfo_2.zip"
echo "Unzip the content until you have: models/M_barkeri_model_iMG746/mba_model.xml"
echo "Finall run: prepare_M_markeri_iMG746_model.py" 
