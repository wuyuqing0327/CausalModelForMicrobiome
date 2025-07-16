# import undetected_chromedriver.v2 as uc
import undetected_chromedriver as uc
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
from selenium.webdriver.common.action_chains import ActionChains

# Initialize undetected ChromeDriver
driver = uc.Chrome()

# Navigate to the Mummichog server
driver.get('https://mummichog-2.appspot.com/')

# Wait for the file upload element to be present
upload_element = WebDriverWait(driver, 50).until(
    EC.presence_of_element_located((By.CSS_SELECTOR, "input[type='file']"))
)
print("Found the file input element.")

upload_element.send_keys(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichog\microbiome_002.txt')
print("File uploaded successfully.")

# Fill out the email field
email_element = WebDriverWait(driver, 50).until(
    EC.presence_of_element_located((By.NAME, "email"))
)
email_element.send_keys("yuqingw1@uchicagomedicine.org")

# Fill out the adduct_type field
adduct_type_element = WebDriverWait(driver, 50).until(
    EC.presence_of_element_located((By.NAME, "adduct_type"))
)
adduct_type_element.send_keys('["M+H[1+]", "M+Na[1+]", "M[1+]"]')

# Fill out the cutoff field
cutoff_element = WebDriverWait(driver, 50).until(
    EC.presence_of_element_located((By.NAME, "cutoff"))
)
cutoff_element.send_keys("0.05")

# Check the agree to privacy policy checkbox
checkbox_element = WebDriverWait(driver, 50).until(
    EC.element_to_be_clickable((By.NAME, "agreement"))
)
checkbox_element.click()

# Submit the form
submit_element = WebDriverWait(driver, 50).until(
    EC.element_to_be_clickable((By.CSS_SELECTOR, "input[type='submit']"))
)

# Use ActionChains to click the submit button
actions = ActionChains(driver)
actions.move_to_element(submit_element).perform()
time.sleep(2)  # Slightly longer delay to mimic more realistic human interaction
submit_element.click()
print("Submitted successfully.")

# Store the current URL
current_url = driver.current_url

# Wait for the URL to change after submission
WebDriverWait(driver, 100).until(EC.url_changes(current_url))
print("Page has redirected successfully.")

# Wait for a specific element on the new page
WebDriverWait(driver, 300).until(
    EC.presence_of_element_located((By.XPATH, "//*[contains(text(),'Results')]"))
)
print("Results page loaded successfully.")

# Close the browser
driver.quit()
