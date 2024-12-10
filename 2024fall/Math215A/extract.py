import os
import sys
from PyPDF2 import PdfReader, PdfWriter
import subprocess

def extract_pages(input_pdf, output_pdf, start_page, end_page):
    try:
        reader = PdfReader(input_pdf)
        writer = PdfWriter()

        for page_num in range(start_page - 1, end_page):
            writer.add_page(reader.pages[page_num])

        with open(output_pdf, "wb") as output_file:
            writer.write(output_file)

        print(f"Pages {start_page}-{end_page} have been extracted to {output_pdf}")
    except Exception as e:
        print(f"Error: {e}")

def check_conda_env_for_pypdf():
    try:
        # Get list of all Conda environments
        env_list_output = subprocess.check_output(["conda", "env", "list"], text=True)
        envs = [line.split()[0] for line in env_list_output.splitlines() if line and not line.startswith("#")]

        print("Checking for `pypdf` in Conda environments...")

        for env in envs:
            try:
                result = subprocess.run(["conda", "list", "pypdf", "-n", env], capture_output=True, text=True)
                if "pypdf" in result.stdout:
                    print(f"`pypdf` found in environment: {env}")
                    return env
            except Exception as e:
                print(f"Could not check environment {env}: {e}")

        print("`pypdf` not found in any Conda environment.")
        print("Available environments:")
        for env in envs:
            print(f"- {env}")
    except Exception as e:
        print(f"Error checking Conda environments: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_pdf_pages.py <input_pdf>")
        sys.exit(1)

    input_pdf = sys.argv[1]
    output_pdf = "extracted_pages.pdf"
    start_page = 265    
    end_page = 281

    print("Step 1: Checking Conda environments for `pypdf`...")
    check_conda_env_for_pypdf()

    print("Step 2: Extracting pages...")
    extract_pages(input_pdf, output_pdf, start_page, end_page)

