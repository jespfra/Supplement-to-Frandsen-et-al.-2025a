print("Hello World")

# import os
# run1 = os.system.subprocess(["ls", "-l"])
# print(run1)


import subprocess

# Function to run Julia script
def run_julia_script(script_path):
    command = f"julia {script_path}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, error = process.communicate()

    # Check if the process exited with an error
    if process.returncode != 0:
        print(f"Error running {script_path}:")
        print(error)
    else:
        print(f"Output from {script_path}:")
        print(output)


# SMB Linear scripts
# run_julia_script("Linear/SMB/LRM/CADET-Julia-LRMLinear-SemiAnalytical.jl")
# run_julia_script("Linear/SMB/LRM/CADET-Julia-LRMLinear.jl")
# run_julia_script("Linear/SMB/LRM/CSS/CADET-Julia-LRMLinear-SemiAnalytical.jl")
# run_julia_script("Linear/SMB/LRM/CSS/CADET-Julia-LRMLinearCSS.jl")


# SMB Langmuir scripts
# run_julia_script("Langmuir/SMB/LRMP/CADET-Julia-LRMPLangmuir.jl")
# run_julia_script("Langmuir/SMB/LRMP/CSS/CADET-Julia-LRMPLangmuir-SemiAnalytical.jl")
# run_julia_script("Langmuir/SMB/LRMP/CSS/CADET-Julia-LRMPLangmuirCSS.jl")


# SMB SMA scripts
# run_julia_script("SMA/SMB/LRMP/CADET-Julia-LRMPSMA.jl")
# run_julia_script("SMA/SMB/LRMP/CSS/CADET-Julia-LRMPSMA-SemiAnalytical.jl")
# run_julia_script("SMA/SMB/LRMP/CSS/CADET-Julia-LRMPSMACSS.jl")
# run_julia_script("SMA/SMB/GRM/CADET-Julia-GRMSMA.jl")
# run_julia_script("SMA/SMB/GRM/CSS/CADET-Julia-GRMSMA-SemiAnalytical.jl")
# run_julia_script("SMA/SMB/GRM/CSS/CADET-Julia-GRMSMACSS.jl")


