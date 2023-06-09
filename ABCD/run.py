import argparse
import os
# import torch
import random
import numpy as np
from ABCD.preprocess.preprocessing import *
from ABCD.train.trainer import *

def seed_everything(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    # cudnn.benchmark = True
    # cudnn.deterministic = True

def parse_args(jupyter=False):
    parser = argparse.ArgumentParser()
    ######################## DATA Loading ########################
    parser.add_argument('--preprocess', type = int, default=1, help='preprocess the data or not')
    parser.add_argument('--group', type=str, default='anxiety', help = 'Which group to include, supported groups include adhd, bipolar, anxiety, panic, ocd')
    parser.add_argument('--label_path', type=str, 
                        default = '/gpfs/slayman/pi/gerstein/jjl86/DATA/ABCD/labels/2022_03_18_label_anxiety_group_versus_nonclinical_controls.csv,
                        help = 'path of the label files, 1 for disease, 0 for ctrl')
    parser.add_argument('--raw_data_path', type=str,
                        default = '/gpfs/slayman/pi/gerstein/jjl86/DATA/ABCD/aurora01_combined',
                        help = 'raw data path')
    parser.add_argument('--cov_path', type=str,
                        default = '/gpfs/slayman/pi/gerstein/jjl86/DATA/ABCD/covariates_processed.csv',
                        help = 'metadata for subjects, storing the demographical data')
    parser.add_argument('--out_path', type=str,
                        default = '/gpfs/slayman/pi/gerstein/jjl86/DATA/ABCD/processed_data')
    parser.add_argument('--save_path', type=str,
                        default = '/gpfs/slayman/pi/gerstein/jjl86/DATA/ABCD/results')
    
    parser.add_argument('--add_genomics', type=int,
                        default=1)
    
    parser.add_argument('--X_path', type=str, help = 'path of the processed X data, only specify this if you have run preprocess function',
                        default= '/home/yl2428/ABCD/processed_data/new_genomic_all_9000_adhd_non_clinical_2_3_X.npy' )
    parser.add_argument('--Y_path', type=str, help = 'path of the processed Y data, only specify this if you have run preprocess function',
                        default= '/home/yl2428/ABCD/processed_data/new_genomic_all_9000_adhd_non_clinical_2_3_Y.npy' )
    parser.add_argument('--subject_path', type=str, help ='path to save the subject of the splits', default = './')

    
    ####################### Preprocessing ########################
    parser.add_argument('--filter_threshold', type=int, default=9000, 
                        help = 'QC control parameter. Normally for a two day frame, 8000 to 10000 will be good choices')
    parser.add_argument('--subseted_days', type=list,
                        default=[1,2],
                        help = 'which days are incoporated into the pipeline')

    ######################## Model ########################
    parser.add_argument('--pretrain', type = int, default = 0)
    parser.add_argument('--freeze_epoches', type=int, default=0)
    parser.add_argument('--model', type= str, default = 'XceptionTimePlus')
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--Grad_CAM_size', type=int, default=32, 
                            help ='number of samples for grad cam to compute average over')
    parser.add_argument('--Grad_CAM_layer', type = int, default = 0,
                            help = 'Grad CAM targeting layer')
    parser.add_argument('--step_importance_steps', type = int, default=60,
                            help = "resolution of step importance, default is 60 measurements (1hr)")
    ###################### Experiments Settings ############
    parser.add_argument('--n_tests', type = int, default=1,
                            help = "number of exp to run)")
    parser.add_argument('--epochs', type = int, default=12,
                            help = "number of epochs to train)")
    parser.add_argument('--test_size', type = float, default=0.2,
                            help = "size_of_the_test")
    parser.add_argument('--class_balance', type = int, default = 0,
                            help = "Perform class balancing or not")
    parser.add_argument('--grad_cam', type = int, default = 0,
                            help = "Perform Grad CAM or not")
    parser.add_argument('--step_importance', type = int, default = 0,
                            help = "Perform step_importanc or not")
    parser.add_argument('--feature_importance', type = int, default = 0,
                            help = "Perform feature_importance or not")
    


    if(jupyter):
        args = parser.parse_args(args = [])
    else:
        args = parser.parse_args()
    
    # args.num_workers = min(args.num_workers,args.batch_size)
    # if os.environ["LOCAL_RANK"] is not None:
    #     args.local_rank = int(os.environ["LOCAL_RANK"])
    config = {}
    for key, value in vars(args).items():
        config[key] = value
    return config

def main():
    config = parse_args()
    if config['preprocess']:
        X_out, Y = process_data(config)
    # else:
    #     X_out = np.load(config['X_path'])
    #     Y = np.load(config['Y_path'])
    #     clf = ABCDTrainer(config)
    #     clf.fit(X_out, Y, n_tests = config['n_tests'], n_epochs = config['epochs'], 
    #     test_size = config['test_size'], grad_cam = config['grad_cam'], step_importance = config['step_importance'],
    #     feature_importance = config['feature_importance'])
        

if __name__ == '__main__':
    main()    