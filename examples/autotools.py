import os, re, shutil, subprocess, glob, sys
import numpy as np

def dist(rec, N, s):
    import numpy as np
    y = np.cos(np.pi * np.arange(N) / (N - 1))
    x = rec - s * y / np.sqrt(1 + ((s / (N - 1))**2) - y**2)
    x = np.round(x, 0)
    x = list(map(str, x.astype(int)))
    print('range distribution around ',rec,x)
    return x

def build_nek(folder_path,usr_file,verbose=False):
    from subprocess import Popen, PIPE, STDOUT
    from pathlib import Path
    source_root = os.environ['NEKSTAB_SOURCE_ROOT']
    cwd = os.getcwd() + '/' + folder_path
    print("Compiling nek5000...")
    print(f'    Using working directory "{cwd}"')
    print(f'    Using .usr file "{usr_file}"')
    my_env = os.environ.copy()
    makenek_in = Path(source_root) / "bin" / "mks"
    logfile = Path(cwd) / "build.log"
    proc = Popen([makenek_in, "clean"], cwd=cwd, env=my_env, stdin=PIPE, text=True)
    proc.communicate(input="Y\n")
    proc.wait()
    proc = Popen([makenek_in, usr_file], cwd=cwd, env=my_env, stdin=PIPE, stderr=STDOUT)
    proc.wait()

    if proc.returncode != 0:
        with open(logfile, "r") as file:
            text = file.read()
        print(text)
        exit(-1)

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.exists(d):
            try:
                shutil.rmtree(d)
            except Exception as e:
                print(e)
                os.unlink(d)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def cSZ(infile, outfile,params):
    with open(infile, "r") as f:
        lines = f.readlines()
    # Substitute all the variables
    for key, value in params.items():
        print("Modf.",key, "to",value,'in',infile)
        if value:
            lines = [
                re.sub(
                    r"(.*\bparameter\b.*\b{0} *= *)\S+?( *[),])".format(key),
                    r"\g<1>{0}\g<2>".format(value),
                    l,
                    flags=re.I,
                )
                for l in lines
            ]
    with open(outfile, "w") as f:
        f.writelines(lines)


    
def check_keyword_in_file(folder_path, file_name, keyword):
    file_path = os.path.join(folder_path, file_name)
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            file_contents = file.read()
            if keyword in file_contents:
                print(f"The keyword '{keyword}' was found in the file '{file_name}'.")
                return True
            else:
                print(f"The keyword '{keyword}' was not found in the file '{file_name}'.")
                return False
    else:
        print(f"The file '{file_name}' does not exist in the folder '{folder_path}'.")
        return False

def write_to_file(filename, content, target_line):
    with open(filename, 'r') as file:
        lines = file.readlines()
    #content_str = ' '.join(map(str, content))
    lines[target_line - 1] = content + '\n'
    with open(filename, 'w') as file:
        file.writelines(lines)

def compile_neks(casename):
    command = ["makeneks", casename]
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    output, error = process.communicate(input=b'y\n')  # Send 'y' and newline ('\n')

    if process.returncode != 0:  # If the command failed...
        print(f"Compilation for {casename} failed.")
        if error:
            print("Error message:", error.decode('utf-8'))
    else:
        print("Command output:", output.decode('utf-8'))

def c_pf(infile, outfile, opts):
    import configparser
    parfile = configparser.ConfigParser()
    parfile.read(infile)
    for section, name_vals in opts.items():
        for name, val in name_vals.items():
            print("In",section,':',name,'set to',val,'in',infile)
            parfile.set(section, name, val)
    with open(outfile, "w") as f:
        parfile.write(f)

def check_last_value(filename, tolerance):
    with open(filename, 'r') as file:
        last_line = file.readlines()[-1]   # Read the last line of the file
    last_value = float(last_line.split()[-1])  # Split the line into parts and convert the last part to float
    return last_value < tolerance

def delete_files(pattern):
    [os.remove(file) for file in glob.glob(pattern)]
    
def copy_bf(file,residu_file,oldbfs,tolerance):    
    if os.path.isfile(file) and os.path.isfile(residu_file):
        print(f"Found '{file}' and '{residu_file}'")
        with open(residu_file, 'r') as f:
            lines = f.readlines()
            last_value = float(lines[-1].split()[-1])
            if last_value < tolerance:
                print(f" Last value of '{residu_file}' is {last_value} < {tolerance}")
                curr_folder = os.path.basename(os.getcwd())
                dest_path = os.path.join(oldbfs, curr_folder)
                print(f" Copying '{os.getcwd()+'/'+file}' to '{dest_path}'")
                shutil.copy(file, dest_path)

def check_job_status(job_name):
    result = subprocess.run(['squeue', '-u', 'rvpo014'], capture_output=True, text=True)
    job_status = job_name in result.stdout
    if job_status:
        print(f'Job {job_name} is running.')
        print('')
        print('')
        sys.exit()
    else:
        print(f'Job {job_name} is not running.')
    return job_status

def resubmit_job(pbs_file,folder_path,job_name):
    print(f" Job name: {job_name}")
    pbs_file_path = os.path.join(folder_path, pbs_file)
    print(f" Adjusting file '{pbs_file_path}'")
    with open(pbs_file_path, 'r') as file:
        filedata = file.read()
    match = re.search(r'(#SBATCH -J )\S+', filedata)
    if match:
        old_line = match.group(0)  # Get the entire matched line
        new_line = f"#SBATCH -J {job_name}"
        filedata = filedata.replace(old_line, new_line)
        print(f"Replacing '{old_line}' with '{new_line}'")
    with open(pbs_file_path, 'w') as file:
        file.write(filedata)
    current_dir = os.getcwd()
    os.chdir(folder_path)
    try:
        result = subprocess.run(["sbatch", pbs_file], check=True, stdout=subprocess.PIPE)
        print(f"Job {job_name} submitted successfully.")
        print("Command output:", result.stdout.decode('utf-8'))  # Print the command output
    except subprocess.CalledProcessError:
        print(f"Job submission for {job_name} failed.")
    os.chdir(current_dir)
    
def get_closest_filename(target_dir, reference):
    file_names = os.listdir(target_dir)
    reference = int(reference)  # assuming reference is a numerical value

    print(f' Looking for closest match to {reference} in {target_dir}')
    print(f' Found {file_names} files in {target_dir}')

    closest_diff = float('inf')
    closest_filename = None

    for file_name in file_names:
        # assuming filenames are numbers
        file_number = int(file_name)
        diff = abs(reference - file_number)

        if diff < closest_diff:
            closest_diff = diff
            closest_filename = file_name

    print(f' Closest match to {reference} is {closest_filename} with a difference of {closest_diff}')
    return closest_filename

def adjust_and_submit_job(folder_path, job_name):
    """Adjust the PBS file and submit the job."""
    pbs_file_path = os.path.join(folder_path, pbs_file)
    with open(pbs_file_path, 'r') as file:
        filedata = file.read()

    filedata = filedata.replace('JOBNAME', job_name)
    filedata = filedata.replace('CASENAME', f'"{cn}"')

    with open(pbs_file_path, 'w') as file:
        file.write(filedata)

    try:
        result = subprocess.run(["sbatch", pbs_file], check=True, stdout=subprocess.PIPE)
        print(f"Job {job_name} submitted successfully.")
        print("Command output:", result.stdout.decode('utf-8'))
    except subprocess.CalledProcessError:
        print(f"Job submission for {job_name} failed.")
        
def adjust_job_2(folder_path,pbs_file):
    """Adjust the PBS file and submit the job."""
    pbs_file_path = os.path.join(folder_path, pbs_file)
    print(pbs_file_path)
    try:
        with open(pbs_file_path, 'r') as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            sline = line.strip()
            if sline == '##SBATCH --qos qos_cpu-dev':
                lines[i] = '#SBATCH --qos qos_cpu-dev\n'
            if sline == '#SBATCH --time=10:00:00':
                lines[i] = '#SBATCH --time=01:00:00\n'
            if sline == '#SBATCH --time=20:00:00':
                lines[i] = '#SBATCH --time=01:00:00\n'

        with open(pbs_file_path, 'w') as file:
            file.writelines(lines)
    
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def check_time(filename, which_time='final'):
    # print(f"Checking {which_time} time in {filename}")
    with open(filename, 'r') as file:
        lines = file.readlines()
        #print(lines)
        if not lines:  # if the file is empty
            raise ValueError(f"The file '{filename}' is empty.")
        num = int(lines[0].split()[0])+1  # get the number on the first line
        #print(f"  Number of lines: {num}")
        #print(f"  Number of lines in file: {len(lines)}")
        if len(lines) > num: 
            if which_time == 'final':
                line = lines[-1]
            elif which_time == 'initial':
                line = lines[num]  # skip the first line + 'num' lines
            value = float(line.split()[0])
        else:
            print(f"  File '{filename}' does not contain enough lines.")
            value = 0.0
    #print(f"  Time: {value}")
    return value

def reset_his_file(filename):
    temp_filename = filename + '.tmp'
    with open(filename, 'r') as file, open(temp_filename, 'w') as temp_file:
        num = int(next(file).split()[0])  # read the first line
        temp_file.write(str(num) + '\n')  # write the first line to the temp file
        for _ in range(num):  # copy 'num' lines
            temp_file.write(next(file))
    os.remove(filename)  # delete the original file
    shutil.move(temp_filename, filename)  # move the temp file to the original file's location

def extract_time_from_binary(filename):
    result = subprocess.run(['head', '-1', filename], stdout=subprocess.PIPE)
    first_line = result.stdout.decode('utf-8', errors='ignore')
    words = first_line.split()
    time = float(words[7])  # adjust this index based on where the value is located
    # print(f"Time: {time}")
    return time

def append_files(pattern, output_file):
    files = sorted(glob.glob(pattern))
    with open(output_file, 'w') as outfile:
        for file in files:
            with open(file, 'r') as infile:
                shutil.copyfileobj(infile, outfile)