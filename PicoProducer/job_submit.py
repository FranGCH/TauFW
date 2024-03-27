import os
import sys
import datetime
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--e ", type=str, dest="era", help="which era")
parser.add_argument("--c ", type=str, dest="channel", help="which channel", default='mutau')
parser.add_argument("--s ", nargs='+', dest="systematics", help="which systematics", default = ['central'])
parser.add_argument("--C ", action='store_true', dest="clean", help="want to clean root files", default=False)
options = parser.parse_args()
era = options.era
channel = options.channel
systematics = options.systematics
clean = options.clean

def run_cmd(which_cmd):

    def run_pico(cmd):
        ret = 0
        if 'hadd' in cmd:
            print('running hadd')
            print(cmd)
            os.system(cmd)
        elif 'status' in cmd:
            if not os.system(f'{cmd} | grep MISS'):
                cmd = cmd.replace('status', 'resubmit')
                print('jobs are missing, resubmitting')
                print(cmd)
                os.system(cmd)
                ret = 1
            elif not os.system(f'{cmd} | grep PEND'):
                print(cmd)
                print('jobs are on pending state')
                ret = 1
            else:
                print('jobs are done')
        else:
            if not os.system(f'{cmd} | grep exist'):
                print('jobs were already submitted, check the status')
                cmd = cmd.replace(which_cmd, 'status')
                print(cmd)
                ret = run_pico(cmd)
            else:
                print(cmd)
                os.system(cmd)
                ret = 1
        return ret
    common_cmd = f'pico.py {which_cmd}  -c {channel} -y {era}  -E jec=False useT1=False'
    final_ret = 0
    for sys in systematics:#['central']:#, '_LTFUp', '_LTFDown', '_JTFUp', '_JTFDown', 'tauES']:
        if 'central' in sys:
            ret = run_pico(common_cmd)
            if ret == 1: final_ret = 1
        elif 'LTF' in sys:
            cmd = common_cmd + f' -s DY TT -t {sys}'
            if 'Up' in sys:
                cmd += ' -E ltf=1.03'
            else:
                assert('Down' in sys)
                cmd += ' -E ltf=0.97'
            ret = run_pico(cmd)
            if ret == 1: final_ret = 1
        elif 'JTF' in sys:
            cmd = common_cmd + f' -s DY TT W*J -t {sys}'
            if 'Up' in sys:
                cmd += ' -E jtf=1.10'
            else:
                assert('Down' in sys)
                cmd += ' -E jtf=0.90'
            ret = run_pico(cmd)
            if ret == 1: final_ret = 1
        elif 'tauES' in sys:
               taues_vars = [0.970, 0.972, 0.974, 0.976, 0.978, 0.980, 0.982, 0.986, 0.988, 0.990, 0.992, 0.994, 0.996, 0.998, 1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.012, 1.014, 1.016, 1.018, 1.020, 1.022,1.024, 1.026, 1.028, 1.030]
               for taues_var in taues_vars:
                   taues_var_with_p = str(taues_var).replace('.', 'p')
                   if len(taues_var_with_p) != 5:
                       taues_var_with_p += '0'
                   cmd = common_cmd + f' tes={taues_var} -s DY -t _TES{taues_var_with_p}'
                   ret = run_pico(cmd)
                   if ret == 1: final_ret = 1
    return final_ret


first_submit = False
while True:
    if not first_submit:
        ret = run_cmd('submit')
        print('return value: ', ret)
        if ret == 0:
            break
        first_submit = True
    time.sleep(600)
    ret = run_cmd('status')
    print('return: ', ret)
    print(datetime.datetime.now())
    if ret == 0:
        break
run_cmd('hadd clean' if clean else 'hadd')
