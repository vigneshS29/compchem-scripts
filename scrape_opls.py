import os
import time
import shutil
import re
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# === Load SMILES from Excel ===
df = pd.read_excel('Electrolytes.xlsx')
smiles_list = df['rdkit_smiles'].dropna().tolist()
print(len(smiles_list), "SMILES found in the file.")

# === Directories ===
download_dir = os.getcwd()
output_dir = os.path.join(download_dir, "output_lmp")
os.makedirs(output_dir, exist_ok=True)

# === Chrome Driver Setup ===
options = webdriver.ChromeOptions()
prefs = {
    "download.default_directory": download_dir,
    "download.prompt_for_download": False,
    "safebrowsing.enabled": True
}
options.add_experimental_option("prefs", prefs)

driver = webdriver.Chrome(
    service=Service(ChromeDriverManager().install()),
    options=options
)

wait = WebDriverWait(driver, 60)

# === Function to sanitize filenames ===
def sanitize_filename(name):
    return re.sub(r'[<>:"/\\|?*]', '_', name)

# === Process each SMILES ===
for smiles in smiles_list:
    try:
        print(f"Processing: {smiles}")
        
        # Step 1: Open LigParGen
        driver.get("https://traken.chem.yale.edu/ligpargen/")
        
        # Step 2: Enter SMILES
        wait.until(EC.presence_of_element_located((By.ID, "smiles"))).clear()
        driver.find_element(By.ID, "smiles").send_keys(smiles)
        
        # Step 3: Click "Submit Molecule"
        driver.find_element(By.XPATH, "//button[text()='Submit Molecule']").click()
        
        # Step 4: Wait for "LAMMPS" button to appear
        lammps_button = wait.until(
            EC.element_to_be_clickable((By.XPATH, "//input[@type='submit' and @value='LAMMPS']"))
        )
        lammps_button.click()
        
        # Step 5: Wait for file download
        timeout = time.time() + 60
        downloaded_file = None
        while time.time() < timeout:
            files = [f for f in os.listdir(download_dir) if f.endswith(".lmp")]
            if files and not any(f.endswith(".crdownload") for f in os.listdir(download_dir)):
                downloaded_file = files[0]
                break
            time.sleep(1)
        
        # Step 6: Move and rename
        if downloaded_file:
            safe_name = sanitize_filename(smiles)
            shutil.move(
                os.path.join(download_dir, downloaded_file),
                os.path.join(output_dir, f"{safe_name}.lmp")
            )
            print(f"✅ Saved: {safe_name}.lmp")
        else:
            print(f"❌ Download failed for {smiles}")
    
    except Exception as e:
        print(f"⚠ Error processing {smiles}: {e}")

# Close browser
driver.quit()
