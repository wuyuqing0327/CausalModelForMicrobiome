from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
from selenium.webdriver.common.action_chains import ActionChains


# Install and setup ChromeDriver using webdriver_manager
service = Service(ChromeDriverManager().install())

# Initialize the WebDriver
driver = webdriver.Chrome(service=service)

# Navigate to the Mummichog server
driver.get('https://mummichog-2.appspot.com/')

# Wait for the file upload element to be present
upload_element = WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.CSS_SELECTOR, "input[type='file']"))
)
print("Found the file input element.")

upload_element.send_keys(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichog\microbiome_002.txt')
print("File uploaded successfully.")

# Fill out the email field
email_element = WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.NAME, "email"))
)
email_element.send_keys("yuqingw1@uchicagomedicine.org")

# Fill out the adduct_type field
adduct_type_element = WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.NAME, "adduct_type"))
)
adduct_type_element.send_keys('["M+H[1+]", "M+Na[1+]", "M[1+]"]')

# Fill out the cutoff field
cutoff_element = WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.NAME, "cutoff"))
)
cutoff_element.send_keys("0.05")

# Submit the form
submit_element = WebDriverWait(driver, 20).until(
    EC.element_to_be_clickable((By.ID, "submit"))
)

actions = ActionChains(driver)
actions.move_to_element(submit_element).perform()
time.sleep(1)
submit_element.click()
print("Submitted successfully.")


# Store the current URL
current_url = driver.current_url
# Wait for the URL to change after submission
WebDriverWait(driver, 100).until(EC.url_changes(current_url))
# Now the new page should be loaded
print("Page has redirected successfully.")

# Wait for a specific element on the new page
WebDriverWait(driver, 100).until(
    EC.presence_of_element_located((By.XPATH, "//*[contains(text(),'Results')]"))  # Adjust the text or element as needed
)
print("Results page loaded successfully.")

WebDriverWait(driver, 300).until(
    EC.presence_of_element_located((By.XPATH, "//*[contains(text(),'Results')]"))
)


# Close the browser
driver.quit()
