import os
import logging 
from multiprocessing import Pool
import subprocess

'''
This script is used to transfer the species peak location to mm10 genome location.
'''

def run_hal_liftover(hal_file, query, peak_bed, target, output_psl):
    """
    调用 halLiftover 工具生成 psl 文件。
    
    :param hal_file: HAL 文件路径
    :param query: 查询基因组名称
    :param peak_bed: 输入的 BED 文件路径
    :param target: 目标基因组名称
    :param output_psl: 输出 PSL 文件路径
    """
    cmd = [
        "halLiftover", 
        "--outPSLWithName", 
        hal_file, 
        query, 
        peak_bed, 
        target, 
        output_psl
    ]
    
    # 执行命令并捕获输出
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # 检查命令是否成功运行
    if result.returncode != 0:
        print(f"Error occurred while running halLiftover: {result.stderr}")
    else:
        print("halLiftover executed successfully.")
        print(result.stdout)
    return None


def psl_filter(psl_file, output_bed):
    
    with open (psl_file)as f1, open(output_bed, 'w') as f2:
        for line in f1:
            cols=line.strip().split()
            if (int(cols[17])-int(cols[16])) < 2*(int(cols[13])-int(cols[12])) and (int(cols[13])-int(cols[12])) >30  and 2*(int(cols[17])-int(cols[16])) > (int(cols[13])-int(cols[12])):
                #print (line.strip())
                f2.write(f"{cols[14]}\t{cols[16]}\t{cols[17]}\t{cols[0]}\n")
    return None 

def pipeline(args):
    '''
    Args:
        file_name: prefix for the following output file
        file_path: path of the 10X bed file
    '''
    file_name, file_path = args
    species = file_name[:-2]
    file_dic = {'CY': 'MFU',
                'JT': 'RSI',
                'T': 'zws',
                'QF': 'CSP',
                'M': 'mm10'}

    hal_dic = {'QF':'/home/rsun@ZHANGroup.local/sly_data/coordinate_transfer/maf_hal_file/mm10_CSP.hal',
                'T': '/home/rsun@ZHANGroup.local/sly_data/coordinate_transfer/maf_hal_file/mm10_zws_clean.hal',
                'CY': '/home/rsun@ZHANGroup.local/sly_data/coordinate_transfer/maf_hal_file/5species.hal',
                'JT': '/home/rsun@ZHANGroup.local/sly_data/coordinate_transfer/maf_hal_file/5species.hal'}
    save_dir = f'peaks_mus/{file_name}'
    os.makedirs(save_dir, exist_ok= True)

    # create logger
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        filename=os.path.join(save_dir,'pipeline.log'),  # 日志文件名
                        filemode='a')  # 如果文件已存在，则附加到文件末尾
    
    logging.info(f'parameter: {file_name}, {file_path}')
    
    
    run_hal_liftover(hal_file = hal_dic[species], 
                     query = file_dic[species], 
                     peak_bed = file_path,
                     target = 'mm10', 
                     output_psl = os.path.join(save_dir, 'tmp.psl'))
    logging.info(f'hal liftover over')
    
    psl_filter(psl_file= os.path.join(save_dir, 'tmp.psl'),
               output_bed = os.path.join(save_dir, f'{file_name}_output.bed'))
    logging.info(f'psl filter over')
    
    return None 
    

def get_parameter():
    file_name_list = os.listdir('/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/peaks')
    args_list  = []
    for file_name in file_name_list:
        if file_name[:-2] == 'M':
            continue
        a = f'/home/rsun@ZHANGroup.local/sly_data/sly_07_exfig/global_grn/peaks/{file_name}'
        b = os.path.join(a,'peaks.bed')
        args_list.append([file_name,b])

    return args_list

if __name__ == '__main__':

    args_list = get_parameter()

    with Pool(processes=40) as pool:
        pool.map(pipeline, args_list)
    print('DONE')