import argparse
import os
# import torch
import random
from ABCD.preprocess.preprocessing import *

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
    parser.add_argument('--group', type=str, default='adhd', help = 'Which group to include, supported groups include adhd, bipolar, anxiety, panic, ocd')
    parser.add_argument('--label_path', type=str, 
                        default = '/home/yl2428/ABCD/labels/2022_03_18_label_adhd_group_versus_nonclinical_controls.csv',
                        help = 'path of the label files, 1 for disease, 0 for ctrl')
    parser.add_argument('--raw_data_path', type=str,
                        default = '/home/yl2428/ABCD/aurora01_combined',
                        help = 'raw data path')
    parser.add_argument('--cov_path', type=str,
                        default = '/home/yl2428/ABCD/covariates_processed.csv',
                        help = 'metadata for subjects, storing the demographical data')
    parser.add_argument('--outpath', type=str,
                        default = '/home/yl2428/ABCD/processed_data')

    
    ####################### Preprocessing ########################
    parser.add_argument('filter_threshold', type=int, default=9000, 
                        help = 'QC control parameter. Normally for a two day frame, 8000 to 10000 will be good choices')
    parser.add_argument('--subseted_days', type=list,
                        default=[1,2],
                        help = 'which days are incoporated into the pipeline')


    ######################## Model ########################
    parser.add_argument('--model', type=str, default="TorchMD_Norm")
    parser.add_argument('--batch_size', type=int, default=32)
    parser.add_argument('--num_interactions', type=int, default=6)
    parser.add_argument('--adaptive_cutoff', action = 'store_true', default=False)
    parser.add_argument('--short_cutoff_lower', type=float, default=0.0)
    parser.add_argument('--short_cutoff_upper', type=float, default=8.0)
    parser.add_argument('--long_cutoff_lower', type=float, default=0.0)
    parser.add_argument('--long_cutoff_upper', type=float, default=10.0)
    parser.add_argument('--otfcutoff', type=float, default=5.0)
    parser.add_argument('--group_center', type=str, default='center_of_mass')
    parser.add_argument('--hidden_channels', type=int, default=128)
    parser.add_argument('--otf_graph', type = bool, default = True, 
                        help = 'on the fly graph construction, only for TorchMDNorm')
    parser.add_argument('--no_broadcast', action='store_true', default=False)
    
    ######################## Optimizer ########################
    parser.add_argument('--learning_rate', type=float, default=1e-3)
    parser.add_argument('--lr_patience', type=int, default=30)
    parser.add_argument('--fp16', default=False, action='store_true')
    parser.add_argument('--gradient_clip', default=False, action='store_true')
    # parser.add_argument('--AMSGrad', default=False, action='store_true')
    parser.add_argument('--warmup_steps', type = int, default = 1000)
    parser.add_argument('--test_interval', type= int, default = 600)
    
    ######################## Training ########################
    parser.add_argument('--inductive', action='store_true', default=False, help='inductive learning, only implemented for HUT Dataset')
    parser.add_argument('--early_stop', action='store_true', default=False, help='early stopping, default patience is 30')
    parser.add_argument('--early_stop_patience', type=int, default=500, help='early stopping, default patience is 30')
    parser.add_argument('--max_epochs', type=int, default=10000)
    parser.add_argument('--rho_tradeoff', type=float, default=.01)
    parser.add_argument('--sample_size', type=int, default=-1)
    parser.add_argument('--train_prop', type=float, default=0.95)
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--regress_forces', type=bool, default=True)
    parser.add_argument('--num_workers', type=int, default=24) ## gpu device count
    parser.add_argument('--local_rank', type=int, default=0) ## gpu device countpu device count
    parser.add_argument('--master_port', type=str, default="1234") ## gpu device countpu device count
    
    ######################## Logging ########################
    parser.add_argument('--wandb', action='store_true', default=False) ## gpu device count
    parser.add_argument('--name', type=str, default=None, help='name of the experiment')
    parser.add_argument('--notes', type=str, default=None, help='description of the experiment')
    parser.add_argument('--tags', type=str, default=None, help='tags of the experiment')
    parser.add_argument('--group', type=str, default=None, help='project name of the experiment')
    parser.add_argument('--api_key', type=str, default=None, help='wandb api key')

    if(jupyter):
        args = parser.parse_args(args = [])
    else:
        args = parser.parse_args()
    
    args.num_workers = min(args.num_workers,args.batch_size)
    # if os.environ["LOCAL_RANK"] is not None:
    #     args.local_rank = int(os.environ["LOCAL_RANK"])
    config = {}
    for key, value in vars(args).items():
        config[key] = value
    return config

def test():
    config = parse_args()
    process_data(config)

if __name__ == '__main__':
    test()    