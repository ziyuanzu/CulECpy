import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import pfba
from cobra.io import write_sbml_model,save_json_model,load_json_model,read_sbml_model
import re
import pandas as pd
import numpy as np
import xlrd
import re
import openpyxl
import sys
sys.path.append(r'./')
import os
from copy import copy, deepcopy
import numpy as np

def exogenous_split(exogenous,meta_key):
    add_model=True
    exogenous_A=[]
    exogenous_B=[]
    print(meta_key)
    for i in exogenous:
        if i.reversibility==True:
            sub_list=i.reaction.split(' <=> ')[0].split(' + ')
        else:
            sub_list=i.reaction.split(' --> ')[0].split(' + ')
        if meta_key in sub_list:
            add_model=False
        print(sub_list)
        if add_model==True:
            exogenous_A.append(i)
        else:
            exogenous_B.append(i)
    print((exogenous_A,exogenous_B))
    return (exogenous_A,exogenous_B)

def add_new_ex(model,meta):
    meta_name=meta.split('_c')[0]
    meta_c=model.metabolites.get_by_id(meta)
    meta_e=Metabolite(meta_name+'_e',
                formula=meta_c.formula,
                name=meta_c.name,
                compartment='e'
                )
    EX_new=Reaction("EX_"+meta_name+"_e")
    EX_new.name=meta_name+' exchange'
    EX_new.subsystem='unknow'
    EX_new.lower_bound=0
    EX_new.upper_bound=1000
    EX_new.add_metabolites({
        meta_e:-1.0
    })
    new_c2e=Reaction(meta_name+"_c2e")
    new_c2e.name=meta_name+'_c2e'
    new_c2e.subsystem='unknow'
    new_c2e.lower_bound=-1000
    new_c2e.upper_bound=1000
    new_c2e.add_metabolites({
        meta_c:-1.0,
        meta_e:1.0
    })
    model.add_reaction(EX_new)
    model.add_reaction(new_c2e)
    return model

def get_new_model(model_name,exogenous,out_files,meta_key,add_product=False,product=None):
    model = cobra.io.load_json_model(model_name)
    #初始化底物摄入
    model.reactions.get_by_id('EX_glc__D_e').bounds=(-10,1000)
    #引入外源反应
    model.add_reactions(exogenous)
    model=add_new_ex(model,meta_key)
    if add_product==True:
        model=add_new_ex(model,product)
    save_json_model(model,out_files+'.json')
    write_sbml_model(model,out_files+'.xml')
    sbml_excel(out_files+'.xml',out_files+'.xlsx')  


def get_new_system(model_name,exogenous,model_A,model_B,meta_key,product):
    (exogenous_A,exogenous_B)=exogenous_split(exogenous,meta_key)
    #A
    get_new_model(model_name,exogenous_A,model_A,meta_key)
    #B
    get_new_model(model_name,exogenous_B,model_B,meta_key,True,product)

def get_new_AB_model(model_origin_name,model_a_name,model_b_name,meta_test,DM_product,model_ab_json):
    model_a=load_json_model(model_a_name+'.json')
    model_b=load_json_model(model_b_name+'.json')
    model_a.reactions.get_by_id('EX_glc__D_e').bounds=(-10,1000)
    model_a.reactions.get_by_id('EX_'+meta_test+'_e').bounds=(0,1000)
    model_b.reactions.get_by_id('EX_'+meta_test+'_e').bounds=(-1000,0)
    model_b.reactions.get_by_id('EX_glc__D_e').bounds=(-10,1000)
    model_a.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds=(0.1,1000)
    model_b.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds=(0.1,1000)
    print("--------A  model----------")
    print(len(model_a.reactions))
    print(len(model_a.metabolites))
    print(len(model_a.genes) )
    print("--------B  model----------")
    print(len(model_b.reactions))
    print(len(model_b.metabolites))
    print(len(model_b.genes) )
    model=deepcopy(model_a)
    for i in model.reactions:
        i.id='A_'+i.id
    for i in model.metabolites:
        i.id='A_'+i.id
        i.name='A_'+i.name
    model_t=deepcopy(model_b)
    # for i in model_t.genes:
    #     i.id=i.id+'_B'
    for i in model_t.reactions:
        i.id='B_'+i.id
    for i in model_t.metabolites:
        i.id='B_'+i.id
        i.name='B_'+i.name
    for i in model_t.reactions:
        model.add_reaction(i)
    e_meta_rn={}
    for meta in model_a.metabolites:
        if meta.compartment=='e' :
            meta_a=model.metabolites.get_by_id('A_'+meta.id)
            e_meta_rn[meta.id+"_1"]=Reaction('A_'+meta.id+'_con_1')
            e_meta_rn[meta.id+"_1"].name = meta.id+' Connection 1'
            e_meta_rn[meta.id+"_1"].subsystem = 'unKnow'
            e_meta_rn[meta.id+"_1"].lower_bound = 0. # This is the default
            e_meta_rn[meta.id+"_1"].upper_bound = 1000. # This is the default
            e_meta_rn[meta.id+"_1"].add_metabolites({
                meta_a:-1.0,
                meta:1.0
                })
            model.add_reactions([e_meta_rn[meta.id+"_1"]])
    for meta in model_b.metabolites:
        if meta.compartment=='e' :
            meta_b=model_t.metabolites.get_by_id('B_'+meta.id)
            e_meta_rn[meta.id+"_2"]=Reaction('B_'+meta.id+'_con_2')
            e_meta_rn[meta.id+"_2"].name = meta.id+' Connection 2'
            e_meta_rn[meta.id+"_2"].subsystem = 'unKnow'
            e_meta_rn[meta.id+"_2"].lower_bound = -1000. # This is the default
            e_meta_rn[meta.id+"_2"].upper_bound = 0. # This is the default
            e_meta_rn[meta.id+"_2"].add_metabolites({
                meta_b:1.0,
                meta:-1.0
                })
            model.add_reactions([e_meta_rn[meta.id+"_2"]])
    origin_model=cobra.io.load_json_model(model_origin_name)
    # origin_model.reactions.get_by_id('EX_xyl__D_e').bounds=(-10.0,1000)
    origin_model.medium
    for i in origin_model.medium:
        model.add_reaction(origin_model.reactions.get_by_id(i))
        ex_meta=i.split('EX_')[1]
        model.reactions.get_by_id('A_'+ex_meta+'_con_1').bounds=(-1000,1000)
        model.reactions.get_by_id('B_'+ex_meta+'_con_2').bounds=(-1000,1000)
        print('A_'+ex_meta+'_con_1')
    model.reactions.get_by_id('A_'+meta_test+'_e'+'_con_1').bounds=(-1000,1000)
    model.reactions.get_by_id('B_'+meta_test+'_e'+'_con_2').bounds=(-1000,1000)
    for i in model.medium:
        h=i.split('_')[0]
        if h=='A' or h=='B':
            model.reactions.get_by_id(i).bounds=(0,0)
    model.medium
    #model.medium
    model.add_reaction(model_b.reactions.get_by_id(DM_product))
    save_json_model(model,model_ab_json)

def sbml_excel(modelname,output):
    from cobra.io.dict import model_to_dict, model_from_dict,metabolite_from_dict,gene_from_dict,reaction_from_dict
    from cobra.core import Model
    from cobra.io import read_sbml_model, write_sbml_model
    import pandas as pd
    model=read_sbml_model(modelname)
    a=model_to_dict(model,sort=False)
    #writer = pd.ExcelWriter(modelname[:-4]+'.xlsx')
    writer = pd.ExcelWriter(output)
    print(output)
    #pd.DataFrame(a['compartments']).to_excel(writer,'Sheet1',index=False)
    pd.DataFrame(a['metabolites']).to_excel(writer,'metabolites',index=False)
    pd.DataFrame(a['genes']).to_excel(writer,'genes',index=False)
    df_r=pd.DataFrame(a['reactions'])
    df_r['reaction_eq'] = df_r.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=False))
    df_r['reaction_eq_name'] = df_r.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=True))
    #del df_r['metabolites'] #生成反应方程列后将代谢物列删除，可保留以方便再次读取
    df_r.to_excel(writer,'reactions',index=False)
    #df_r.to_excel(writer,'reactions',columns=["id","name","reaction","lower_bound","upper_bound","gene_reaction_rule"],index=False)
    writer.save()

def pathway (model,fluxes,outputfile_name):
    import csv
    import pandas as pd
    flux = open(outputfile_name, 'w')
    for r,v in fluxes.iteritems():
        if abs(v)>1e-6:
            flux.write(r +'\t' + str(round(v,4)) + '\t' + model.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True) +'\n')
    flux.close()

def add_new_ex_special(model,meta_name,meta_c,meta_e):
    meta_c=model.metabolites.get_by_id(meta_c)
    meta_e=Metabolite(meta_e,
                formula=meta_c.formula,
                name=meta_c.name,
                compartment='e'
                )  
    #meta_e=model.metabolites.get_by_id(meta_e)
    EX_new=Reaction("EX_"+meta_name+"_e")
    EX_new.name=meta_name+' exchange'
    EX_new.subsystem='unknow'
    EX_new.lower_bound=-1000
    EX_new.upper_bound=1000
    EX_new.add_metabolites({
        meta_e:-1.0
    })
    new_c2e=Reaction(meta_name+"_c2e")
    new_c2e.name=meta_name+'_c2e'
    new_c2e.subsystem='unknow'
    new_c2e.lower_bound=-1000
    new_c2e.upper_bound=1000
    new_c2e.add_metabolites({
        meta_c:-1.0,
        meta_e:1.0
    })
    model.add_reaction(EX_new)
    model.add_reaction(new_c2e)
    return model