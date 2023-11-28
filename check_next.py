#!/usr/bin/env python3
import sys
sys.path.append("../")
from cube import *

if __name__ == "__main__":
    root = os.getcwd()

    try:
        with open(file_ct, 'r') as file:
            ct = float(file.read().strip())
    except FileNotFoundError:
        ct = 2.0
    except Exception as e:
        print(f"An error occurred while reading the file: {str(e)}")
        ct = 2.0

    print('')
    print('Check next:',ct)
    
    bf_exists = os.path.isfile(bf)
    fnewt_exist = os.path.isfile(fnewt)
    spec_d_exists = os.path.isfile(spec_d)
    spec_d_lines = sum(1 for line in open(spec_d)) if spec_d_exists else 0
    spec_a_exists = os.path.isfile(spec_a)
    spec_a_lines = sum(1 for line in open(spec_a)) if spec_a_exists else 0
    budgetf = os.path.isfile(budget)
        
    folder_path='.'
    pf = glob.glob(os.path.join(folder_path, '*.par'))[0]
    print(f"Parameter file: {pf}")
    cn = pf.split('/')[-1].split('.')[0]
    print(f"Case name: {cn}")
    Re = os.path.basename(os.getcwd())
    print(f"Reynolds number: {Re}")

    if ct >= 2 and bf_exists and fnewt_exist and not spec_d_exists:

        if not check_last_value(fnewt,tol):
            raise Exception("BF not converged to tolerance. CHECK!")
        else:
            print(f" Reading file {fnewt}.")
            print(f" Base flow converged to {tol}.")

        job_name = f"{cns}{Re}d"
        check_job_status(job_name)
        print(f"Submitting direct computation.")
        copy_bf(bf,residu_file,oldbfs,tol)
        delete_files('nwt*')
        shutil.copy(pf, pf+'_b')
        shutil.copy('logfile', 'logfile'+'_b')
        c_pf(pf,pf, {'GENERAL':{'userParam01':'3.1'}})
        resubmit_job(pbs_file,folder_path,job_name)
        
    elif ct >= 2.1 and bf_exists and fnewt_exist and spec_d_lines >= cvgd and not budgetf:
   
        print(f" In ct >= 2.1 running PKE analysis.")
        delete_files('KRY*')
        job_name = f"{cns}{Re}pk"
        check_job_status(job_name)
        print(f"Direct finished, submitting PKE...")
        shutil.copy(pf, pf+'_pk')
        shutil.copy('logfile', 'logfile'+'_pk')
        c_pf(pf,pf, {'GENERAL':{'userParam01':'4.1'}})
        adjust_job_2(folder_path,pbs_file)
        print(f"Recompiling nek5000:")
        build_nek(folder_path, usr_file=cn)
        resubmit_job(pbs_file,folder_path,job_name)
    
    elif ct >= 3 and bf_exists and fnewt_exist and spec_d_exists and not spec_a_exists:

        job_name = f"{cns}{Re}a"
        check_job_status(job_name)
        print(f"Direct computation finished, submitting adjoint.")
        shutil.copy(pf, pf+'_d')
        shutil.copy('logfile', 'logfile'+'_d')
        c_pf(pf,pf, {'GENERAL':{'userParam01':'3.2'}})
        resubmit_job(pbs_file,folder_path,job_name)
        
    elif ct >= 4 and bf_exists and fnewt_exist and spec_d_lines >= cvgd and spec_a_lines >= cvgd and not budgetf:
    
        delete_files('KRY*')
        job_name = f"{cns}{Re}p"
        check_job_status(job_name)
        print(f"Direct and adjoint finished, submitting final.")
        shutil.copy(pf, pf+'_a')
        shutil.copy('logfile', 'logfile'+'_a')
        c_pf(pf,pf, {'GENERAL':{'userParam01':'4'}})
        resubmit_job(pbs_file,folder_path,job_name)
           
    else:
        
        print("No action taken.")
