#! /usr/bin/env python
# -*- coding: utf-8 -*-
import chimera

# use 'rc' as shorthand for runCommand (runCommand is the module which print command to chimera)
from chimera import runCommand as rc 
# for emitting status messages
from chimera import replyobj 

data_dir = "/Users/rechique/Documents/IGEM/Dry Lab/Prototype II/Structure modelling"


"""
PAM_out это когда PAM сайты двух касов находятся на разных цепях ДНК и расположены снаружи от Касов
фотография для понимания: http://parts.igem.org/wiki/images/0/06/Peking-44.png
В стрктуре Каса приняты следующие аннотации:
Chain A - белок
Chain B - гРНК
Chain C - одна цепь ДНК
Chain D - вторая цепь ДНК
"""
    #center_left = 11
    #center_right = 12
def PAM_out_SpCas9_complex(lenght, height): 
    #length - distance between two PAM-sites
    #site - site of connecteion stem to our PDB Cas DNA, it is from -18 to ...
    #height is lenght of space between site and site of connection of the second cas
    """
    тут происходит первичная подготовка модели к дальнейшим операциям. В частности, нуклеотиды ДНК переиндексируются так, что первый нуклеотида PAM имел индекс 0. Нуклеотид слева от него будет -1, а справа +1)
    length - расстояние между PAM-сайтом(точнее комплиментарной ему областью) и местом прикрепления 
    """
    site = -lenght + 3
    height = height - 2
    
    DNA_Cas_complex_length = 41  #length of DNA in model of Cas
    
    #future new start index for dna chain C
    Cas_1_chain_C_end = -18 
    #future new start index for dna chain D
    Cas_1_chain_D_start = -22 
    
    #Download needed structures 
    Cas1 = chimera.openModels.open('%s/5Y36.pdb' % (data_dir), type='PDB')
    #download stem-loop strcture. length of this sequence is 111-(-17)+1=128   
    # -17+x=111 => x=111+17=128; the first index is -17, the second is 128
    Stem = chimera.openModels.open('%s/kirill_stem_without_H.pdb' % data_dir, type='PDB') 
    Cas2 = chimera.openModels.open('%s/5Y36.pdb' % (data_dir), type='PDB')
    
    
    rc('changechains A G #1')  #rename chain A of stem to chain G for comfortable manipulation 
    
    rc('resrenumber %d #0:.C' % (Cas_1_chain_C_end)) 
       
    # delete one chain of Cas (chain D)
    rc('delete #0:.D')   

    # индексируем наш стем-луп так, чтобы при склейке с ДНК Каса индексация сохранилась 
    rc('resrenumber %d #1:.G' % site) 
    # rc('delete #1:%d.G@HO5\'' % site)
    
    # connect dcas9 and ssdna
    rc('match #0:%d.C,%d.C@C5\',C4\',C3\',C2\',C1\' #1:%d.G,%d.G@C5\',C4\',C3\',C2\',C1\'' % (site, site+1, site+1, site))
    
    for i in range(-(Cas_1_chain_C_end - site)):
        rc('delete #0:%d.C' % (Cas_1_chain_C_end+i))
    rc('delete #1:%d.G' % site)
    
    rc('combine #0,1')
    rc('delete #0,1')
    # rc('bond #3:%d.C@C5\' #3:%d.G@OP2' % (site, site+1))
  
    
    rc('changechains A,B a,b #2')
    rc('changechains C,D c,d #2')
    
    '''
    height - длинна комплементарной области стебелька до места присоеденения второго Cas-белка
    base_of_stem1 - индекс первого нуклеотида стебелька
    base_of_stem2 – индекс нуклеотида комплиментарного base_of_stem1
    В случае шпильки Кирилла комплементарная область начинается на ПЯТОМ нуклеотиде с начала (т.е. base_of_stem1=site+4,base_of_stem2=site+125).
    '''
    
    base_of_stem1 = site + 4
    base_of_stem2 = site + 125
    
    c_renum = base_of_stem1+height-17 #это расчет переменной для переиндексации цепей с и d. Для цепи с индекс считается так, чтобы три первых нукелотида цепи со стороны склеивания пересекались по индексам с шпилькой.
    d_renum = base_of_stem2-height-23  #Цепь d индексируется с другого конца, поэтому формула для расчета немного другая. Однако итог похожий - со стороны склеивания три нуклеотида по индекса пересекаются с шпилькой.
     
    c_mat = base_of_stem1 + height
    n_mat = base_of_stem2 - height

    # -13
    # 81
    rc('resrenumber %d #2:.c' % c_renum)
    rc('resrenumber %d #2:.d' % d_renum)
    rc('delete #3:%d-%d.G' % (base_of_stem1+height, base_of_stem2-height))
    
    rc('match #2:%d.c,%d.c,%d.c,%d.d,%d.d,%d.d@C5\',C4\',C3\',C2\',C1\' #3:%d.G,%d.G,%d.G,%d.G,%d.G,%d.G@C5\',C4\',C3\',C2\',C1\'' % (c_mat-1, c_mat-2, c_mat-3, n_mat+1, n_mat+2, n_mat+3, c_mat-1, c_mat-2, c_mat-3, n_mat+1, n_mat+2, n_mat+3))
    
    for i in range(1, 17):
        rc('delete #2:%d.c' % (c_renum+i))
        rc('delete #2:%d.d' % (d_renum+42-i))
       #  rc('delete #3:%d.G' % (d_renum+42-i))       
   
    rc('combine #2,3') #merge two models in one
    rc('delete #2,3')
       
    # rc('bond #0:%d.c@OP2 #0:%d.G@C5\'' % (c_mat, c_mat-1))
    # rc('bond #0:%d.d@C5\' #0:%d.G@OP2'  % (n_mat, n_mat+1))
    
    # rc('minimize nsteps 10')
    rc('write format pdb #0 model_%d_%d.pdb' % (lenght, height)) #save 
    
    rc('close #0')                                                          

for i in range(1, 22):
    PAM_out_SpCas9_complex(i, 9) 
#PAM_out_SpCas9_complex(-17, 7) 
#	rc('write format pdb #3 model_%d_%d.pdb'%(-17, j)) #save   
