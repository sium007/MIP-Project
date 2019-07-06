from Bio import SeqIO
import re
import pandas as pd
from collections import OrderedDict
def get_npa_positions(sequence, pattern):
    sequence = str(sequence.seq)
    positions = re.finditer(pattern, sequence)
    return [p.span()[0] for p in positions]

def get_forger_positions(sequence, valid_positions, pattern_len, f1d=37):
    first_npa = valid_positions[0]
    second_npa = valid_positions[1]
    sid = sequence.id
    sequence = str(sequence.seq)
    
    forger_positions = OrderedDict([
        ('R1' , sequence[first_npa-20]),
        ('R2' , sequence[second_npa-12]),
        ('R3' , sequence[second_npa-3]),
        ('R4' , sequence[second_npa+pattern_len]),
        ('F1' , sequence[first_npa+pattern_len+f1d]),
        ('F2' , sequence[second_npa+pattern_len+1]),
        ('F3' , sequence[second_npa+pattern_len+5]),
        ('SubClass', sid.split(';')[0][2:5])
    ]
    )
    print(sid)
    return forger_positions




if __name__=="__main__":
    excluded=0
    data=[]
    for sequence in SeqIO.parse('newextra.fasta','fasta'):
        
        npa_positions = get_npa_positions(sequence, 'NPA')
        np_positions = get_npa_positions(sequence, 'NP')
        
        found=0
        if len(npa_positions)==2:
            forger_positions = get_forger_positions(sequence, npa_positions,3)
            found=1
        elif len(npa_positions)==3:
            if abs(0-np_positions[0])< 30:
                del npa_positions[0]
            elif abs(len(sequence.seq)-npa_positions[2])<30:
                del npa_positions[2]
            else:
                del npa_positions[1]
            forger_positions = get_forger_positions(sequence,npa_positions,3)
            found=1
        else:
            if len(np_positions)==2:
                second_np_ = str(sequence.seq)[np_positions[1]:np_positions[1]+3]
                print(second_np_=="NPA")
                if second_np_=="NPA":
                    first_np_ = str(sequence.seq)[np_positions[0]:np_positions[0]+3]
                    if first_np_=="NPT" or first_np_=="NPL":
                        forger_positions = get_forger_positions(sequence, np_positions,3,f1d=39)
                        found=1
                    else:
                        forger_positions = get_forger_positions(sequence, np_positions,3)
                        found=1
                else:
                    forger_positions = get_forger_positions(sequence, np_positions,3, f1d=37)
                    found=1
            elif len(np_positions)==3:
                first_np_ = str(sequence.seq)[np_positions[0]:np_positions[0]+3]
                second_np_ = str(sequence.seq)[np_positions[1]:np_positions[1]+3]
                third_np_ = str(sequence.seq)[np_positions[2]:np_positions[2]+3]
                print(first_np_, second_np_, third_np_)
                if first_np_=="NPA":
                    del np_positions[1]
                    forger_positions = get_forger_positions(sequence, np_positions,3)
                    found=1
                elif second_np_=="NPA":
                    dist_1_2= np_positions[1]-np_positions[0]
                    dist_2_3= np_positions[2]-np_positions[1]
                    if dist_1_2>dist_2_3:
                        del np_positions[2]
                        forger_positions = get_forger_positions(sequence, np_positions,3)
                        found=1
                    else:
                        del np_positions[0]
                        forger_positions = get_forger_positions(sequence, np_positions,3)
                        found=1
                elif third_np_=="NPA":
                    dist_2_3 = np_positions[2]-np_positions[1]
                    print('third np',dist_2_3, sequence.id)
                    if dist_2_3>90:
                        del np_positions[0]
                        forger_positions = get_forger_positions(sequence, np_positions,3)
                        found=1
                    else:
                        del np_positions[1]
                        forger_positions = get_forger_positions(sequence, np_positions,3,f1d=39)
                        found=1
                else:
                    del np_positions[1]
                    forger_positions = get_forger_positions(sequence, np_positions,3)
                    found=1
            elif len(np_positions)==1:
                try:
                    spa_npa = [str(sequence.seq).index('SPA'), np_positions[0]]
                    spa_npa.sort()
                    forger_positions = get_forger_positions(sequence, spa_npa,3)
                    found=1
                except:
                    print("Not acquaporin")

            else:
                print("Not acquaporin")
        if found==1:
            # print(forger_positions)
            # print(forger_positions.values())
            data.append(forger_positions)
    
    df = pd.DataFrame(data)
    print(df.columns)

    #print(df.head())
    df.to_csv('newextra.csv', index=False)


##                
##


    