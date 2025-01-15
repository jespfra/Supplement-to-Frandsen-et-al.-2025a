print('Hello World')


import subprocess
import os

# Function to run Julia script
def run_python_script(script_path):
    command = ["python", os.path.abspath(script_path)]
    path = os.path.dirname(os.path.abspath(script_path))
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,cwd=path)
    output, error = process.communicate()

    # Check if the process exited with an error
    if process.returncode != 0:
        print(f"Error running {script_path}:")
        print(error)
    else:
        print(f"Output from {script_path}:")
        print(output)


# SMB Linear scripts
# run_python_script("Linear/SMB/LRM/CADET-FV-LRMLinear.py")
# run_python_script("Linear/SMB/LRM/CADET-DG-LRMLinear.py") 
# run_python_script("Linear/SMB/LRM/CSS/CADET-DG-LRMLinearCSS.py") 
# run_python_script("Linear/SMB/LRM/CSS/CADET-FV-LRMLinearCSS.py") 


# SMB Langmuir scripts
# run_python_script("Langmuir/SMB/LRMP/CADET-FV-LRMPLangmuir.py")
# run_python_script("Langmuir/SMB/LRMP/CADET-DG-LRMPLangmuir.py")
# run_python_script("Langmuir/SMB/LRMP/CSS/CADET-DG-LRMPLangmuirCSS.py")
# run_python_script("Langmuir/SMB/LRMP/CSS/CADET-FV-LRMPLangmuirCSS.py")

# SMB SMA scripts
# run_python_script("SMA/SMB/LRMP/CADET-FV-LRMPSMA.py")
# run_python_script("SMA/SMB/LRMP/CADET-DG-LRMPSMA.py")
# run_python_script("SMA/SMB/LRMP/CSS/CADET-DG-LRMPSMACSS.py")
# run_python_script("SMA/SMB/LRMP/CSS/CADET-FV-LRMPSMACSS.py")

# run_python_script("SMA/SMB/GRM/CADET-FV-GRMSMA.py")
# run_python_script("SMA/SMB/GRM/CADET-DG-GRMSMA.py")
# run_python_script("SMA/SMB/GRM/CSS/CADET-DG-GRMSMACSS.py")
# run_python_script("SMA/SMB/GRM/CSS/CADET-FV-GRMSMACSS.py")







