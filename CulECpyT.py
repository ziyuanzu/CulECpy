import cobra
import numpy as np
import pandas as pd
import pyomo.environ as pyo
from pyomo.environ import *
from pyomo.opt import SolverFactory
import re
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import pfba
import pandas as pd
import numpy as np
import xlrd
import re
import sys
sys.path.append(r'./')
import os
from copy import copy, deepcopy
import numpy as np
from cobra.io import write_sbml_model,load_json_model,save_json_model
import os
from copy import copy, deepcopy
import numpy as np
from cobra.io import write_sbml_model

def new_strip(str):
    while(1):
        if len(str)!=0 and str[0]=='('and str[-1]==')':
            str=str[1:-1]
        else:
            break
    return str

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

def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = new_strip(r.gene_reaction_rule)
    return model

#Building

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

def unpdate_3_kcat_mw(reaction_list,kcat_mw,new_file):#三菌体系更新网络模型酶动力学信息
    reaction_up_kact_MW=pd.DataFrame()
    for i in reaction_list:
        if i.startswith('A_'):
            id=i[2:]
            if id in kcat_mw.index:
                reaction_up_kact_MW.at[i,'kcat_MW']=kcat_mw.at[id,'kcat_MW']
        if i.startswith('B_'):
            id=i[2:]
            if id in kcat_mw.index:
                reaction_up_kact_MW.at[i,'kcat_MW']=kcat_mw.at[id,'kcat_MW']
        if i.startswith('C_'):
            id=i[2:]
            if id in kcat_mw.index:
                reaction_up_kact_MW.at[i,'kcat_MW']=kcat_mw.at[id,'kcat_MW']
    reaction_up_kact_MW.to_csv(new_file)
    return reaction_up_kact_MW

def Model_Solve(model,solver):
    opt = pyo.SolverFactory(solver)
    opt.solve(model)
    return model

def Get_Model_3_Data(model):#获取三菌模型数据
    """Returns reaction_list,metabolite_list,lb_list,ub_list,coef_matrix from model.
    
    Notes: 
    ----------
    *model： is in SBML format (.xml).
    """
    reaction_list=[]
    reaction_list_A=[]
    reaction_list_B=[]
    reaction_list_C=[]
    metabolite_list=[]
    lb_list={}
    ub_list={}
    coef_matrix={}
    for rea in model.reactions:
        reaction_list.append(rea.id)
        if rea.id.startswith('A_'):
            reaction_list_A.append(rea.id)
        if rea.id.startswith('B_'):
            reaction_list_B.append(rea.id)
        if rea.id.startswith('C_'):
            reaction_list_C.append(rea.id)     
        lb_list[rea.id]=rea.lower_bound
        ub_list[rea.id]=rea.upper_bound
        for met in model.metabolites:
            metabolite_list.append(met.id)
            try:
                rea.get_coefficient(met.id)  
            except:
                pass
            else:
                coef_matrix[met.id,rea.id]=rea.get_coefficient(met.id)
    reaction_list=list(set(reaction_list))
    metabolite_list=list(set(metabolite_list))
    return(reaction_list,reaction_list_A,reaction_list_B,reaction_list_C,metabolite_list,lb_list,ub_list,coef_matrix)

def Get_Concretemodel_Need_Data_3_json(model_file,reaction_kcat_MW_file,target=None):#获取三菌模型转化为pyomo数学矩阵形式
    Concretemodel_Need_Data={}
    # reaction_g0=pd.read_csv(reaction_g0_file,index_col=0,sep='\t')
    # Concretemodel_Need_Data['reaction_g0']=reaction_g0
    # metabolites_lnC = pd.read_csv(metabolites_lnC_file, index_col=0, sep='\t')
    # Concretemodel_Need_Data['metabolites_lnC']=metabolites_lnC
    model=cobra.io.load_json_model(model_file)
    try:
        product = model.metabolites.get_by_id(target)
        model.add_boundary(product, type='demand')  # add demand reaction as the objective
    except:
        pass
    cobra.manipulation.modify.convert_to_irreversible(model)
    isoenzyme_split(model)
    reaction_kcat_MW=pd.read_csv(reaction_kcat_MW_file,index_col=0)
    Concretemodel_Need_Data['model']=model
    Concretemodel_Need_Data['reaction_kcat_MW']=reaction_kcat_MW
    [reaction_list,reaction_list_A,reaction_list_B,reaction_list_C,metabolite_list,lb_list,ub_list,coef_matrix]=Get_Model_3_Data(model)
    Concretemodel_Need_Data['reaction_list']=reaction_list
    Concretemodel_Need_Data['reaction_list_A']=reaction_list_A
    Concretemodel_Need_Data['reaction_list_B']=reaction_list_B
    Concretemodel_Need_Data['reaction_list_C']=reaction_list_C
    Concretemodel_Need_Data['metabolite_list']=metabolite_list
    Concretemodel_Need_Data['lb_list']=lb_list
    Concretemodel_Need_Data['ub_list']=ub_list
    Concretemodel_Need_Data['coef_matrix']=coef_matrix
    return (Concretemodel_Need_Data)

def Template_Concretemodel_3(reaction_list=None,metabolite_list=None,coef_matrix=None,reaction_kcat_MW=None,lb_list=None,ub_list=None,reaction_list_A=None,reaction_list_B=None,\
    reaction_list_C=None,set_substrate_ini=False,substrate_name=None,substrate_value=None,substrates=None,substrates_bool=False,\
    set_biomass_ini=False,biomass_value=None,biomass_list=None,biomass_id=None,set_A_biomass_ini=None,set_B_biomass_ini=None,set_C_biomass_ini=None,\
    set_bound=False,set_stoi_matrix=False,set_enzyme_constraint=False,set_part_enzyme_constraint=False,set_obj_single_E_value=False,E_total=None,\
    obj_name=None,obj_target=None,set_obj_value=False,set_target_ini=False,target_value=None,\
    set_obj_V_value=False,set_new_stoi_matrix=False,set_obj_E_value=False,\
    set_ratio=False,set_ratio_bound=False,set_pass_ini=False,set_pass_value=None,\
    set_metabolite=False,set_obj_Met_value=False,metabolites_lnC=None,\
    set_sub_ratio=False,sub_a=None,sub_b=None,sub_ratio=None):
    

    Concretemodel = ConcreteModel()
    Concretemodel.reaction = pyo.Var(reaction_list,  within=NonNegativeReals)
    Concretemodel.z = pyo.Var(reaction_list,  within=pyo.Binary)
    Concretemodel.x = pyo.Var()
    Concretemodel.y = pyo.Var()

    if set_obj_V_value:             
        def set_obj_V_value(m):
            return sum(m.reaction[j] for j in reaction_list)
        Concretemodel.obj = Objective(rule=set_obj_V_value, sense=minimize)  
    
    #Set the maximum flux as the object function
    if set_obj_value:   
        def set_obj_value(m):
            return m.reaction[obj_name]
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_value, sense=maximize)
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_value, sense=minimize)

    #Set the minimum enzyme cost of a pathway as the object function
    if set_obj_E_value:             
        def set_obj_E_value(m):
            return sum(m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index)
        Concretemodel.obj = Objective(rule=set_obj_E_value, sense=minimize)  

    #To calculate the variability of enzyme usage of single reaction.
    if set_obj_single_E_value:             
        def set_obj_single_E_value(m):
            return m.reaction[obj_name]/(reaction_kcat_MW.loc[obj_name,'kcat_MW'])
        if obj_target=='maximize':
            Concretemodel.obj = Objective(rule=set_obj_single_E_value, sense=maximize) 
        elif obj_target=='minimize':
            Concretemodel.obj = Objective(rule=set_obj_single_E_value, sense=minimize)

    if set_sub_ratio:
        def set_sub_ratio(m):
            return m.reaction[sub_a]==4*m.reaction[sub_b]
        Concretemodel.set_sub_ratio = Constraint(rule=set_sub_ratio)

    if set_stoi_matrix:
        def set_stoi_matrix(m,i):
            return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
        Concretemodel.set_stoi_matrix = Constraint( metabolite_list,rule=set_stoi_matrix)

    #Adding flux balance constraints （newFBA）
    # if set_new_stoi_matrix:
    #     def set_new_stoi_matrix(m,i):
    #         for j in reaction_list:
    #             if (i,j) in coef_matrix.keys():
    #                 if j.startswith('A_') and i.startswith('A_') and '_con_' in j:
    #                     coef_matrix[i,j]=(1/m.x)*coef_matrix[i,j]
    #                     # if 'akg_e_con_' in j:
    #                             # print ('A:'+str(coef_matrix[i,j]))
    #                 if j.startswith('B_') and i.startswith('B_') and '_con_' in j:
    #                     coef_matrix[i,j]=(1/m.y)*coef_matrix[i,j]
    #                 if j.startswith('C_') and i.startswith('C_') and '_con_' in j:
    #                     coef_matrix[i,j]=(1/m.z)*coef_matrix[i,j]
    #         return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
    #     Concretemodel.set_new_stoi_matrix = Constraint( metabolite_list,rule=set_new_stoi_matrix)

    #Adding the upper and lower bound constraints of reaction flux
    if set_bound:
        def set_bound(m,j):
            return inequality(lb_list[j],m.reaction[j],ub_list[j])
        Concretemodel.set_bound = Constraint(reaction_list,rule=set_bound) 

    #Set the upper bound for substrate input reaction flux
    if set_substrate_ini and substrates_bool==False:
        def set_substrate_ini(m): 
            return m.reaction[substrate_name] <= substrate_value
        Concretemodel.set_substrate_ini = Constraint(rule=set_substrate_ini)  
    
     #Set the upper bound for substrate input reaction flux
    if set_substrate_ini and substrates_bool==True:
        def set_substrates_ini(m,j): 
            return m.reaction[j]<= substrates[j]
        Concretemodel.set_substrates_ini = Constraint(substrates.keys(),rule=set_substrates_ini)  
    
    if set_biomass_ini:
        def set_biomass_ini(m): 
            return m.reaction[biomass_id] >=biomass_value
        Concretemodel.set_biomass_ini = Constraint(rule=set_biomass_ini)  
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_A_biomass_ini:
        def set_A_biomass_ini(m): 
            return m.reaction['A_'+biomass_id] >=biomass_list['A']
        Concretemodel.set_A_biomass_ini = Constraint(rule=set_A_biomass_ini)  
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_B_biomass_ini:
        def set_B_biomass_ini(m): 
            return m.reaction['B_'+biomass_id] >=biomass_list['B']
        Concretemodel.set_B_biomass_ini = Constraint(rule=set_B_biomass_ini)  
    
    #Set the lower bound for biomass synthesis reaction flux
    if set_C_biomass_ini:
        def set_C_biomass_ini(m): 
            return m.reaction['C_'+biomass_id] >=biomass_list['C']
        Concretemodel.set_C_biomass_ini = Constraint(rule=set_C_biomass_ini)  
    #Set the lower bound for target synthesis reaction flux
    if set_target_ini:
        def set_target_ini(m): 
            return m.reaction[obj_name] ==target_value
        Concretemodel.set_target_ini = Constraint(rule=set_target_ini)  

    #Adding sinngle enzymamic constraints
    if set_enzyme_constraint:
        def set_enzyme_constraint(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list)<= E_total
        Concretemodel.set_enzyme_constraint = Constraint(rule=set_enzyme_constraint)    
    
    #Adding Part A enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_A(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_A)<= E_total
        Concretemodel.set_enzyme_constraint_A = Constraint(rule=set_enzyme_constraint_A)

    #Adding Part B enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_B(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_B)<= E_total
        Concretemodel.set_enzyme_constraint_B = Constraint(rule=set_enzyme_constraint_B)

    #Adding Part C enzymamic constraints 
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_C(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_C)<= E_total
        Concretemodel.set_enzyme_constraint_C = Constraint(rule=set_enzyme_constraint_C)

    #Adding pass meta flux constraints
    if set_pass_ini:
        def set_pass_ini(m,j):
            return m.reaction[j] >= set_pass_value[j]
        Concretemodel.set_pass_ini = Constraint(set_pass_value.keys(),rule=set_pass_ini)  

    return Concretemodel


def thread_FBA(Concretemodel_Need_Data,a,triple,mid_a,mid_b,mid_c,bio_a,bio_b,bio_c,product):#三菌模型进行全局代谢网络模型分析
    #EcoECM=EcoECM(Concretemodel_Need_Data,'DM_B_C06104','maximize',substrates,E_total)\
    #substrates={'EX_glc__D_e_reverse':10}
    x,y,z=triple
    substrates={'EX_glc__D_e_reverse':a}
    substrate_value=10
    substrate_name='EX_glc__D_e_reverse'
    E_total=0.227
    sub_a='EX_glc__D_e_reverse'
    reaction_list_A=deepcopy(Concretemodel_Need_Data['reaction_list_A'])
    reaction_list_B=deepcopy(Concretemodel_Need_Data['reaction_list_B'])
    reaction_list_C=deepcopy(Concretemodel_Need_Data['reaction_list_C'])
    reaction_list=deepcopy(Concretemodel_Need_Data['reaction_list'])
    metabolite_list=deepcopy(Concretemodel_Need_Data['metabolite_list'])
    coef_matrix=deepcopy(Concretemodel_Need_Data['coef_matrix'])
    for i in metabolite_list:
        for j in reaction_list:
                if (i,j) in coef_matrix.keys():
                    if j.startswith('A_') and i.startswith('A_') and '_con_' in j:
                        coef_matrix[i,j]=(1/x)*coef_matrix[i,j]
                        # if 'akg_e_con_' in j:
                                # print ('A:'+str(coef_matrix[i,j]))
                    if j.startswith('B_') and i.startswith('B_') and '_con_' in j:
                        coef_matrix[i,j]=(1/y)*coef_matrix[i,j]
                    if j.startswith('C_') and i.startswith('C_') and '_con_' in j:
                        coef_matrix[i,j]=(1/z)*coef_matrix[i,j]
                        # if 'akg_e_con_' in j:
                #                 print ('B:'+str(coef_matrix[i,j]))
                #     print ('ij:'+str(coef_matrix[i,j]))
    reaction_kcat_MW=deepcopy(Concretemodel_Need_Data['reaction_kcat_MW'])
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    if len(substrates.keys())==1:
            for i in substrates.keys():
                ub_list[i]=substrates[i]
                substrate_name=i
                substrate_value=substrates[i]
            substrates_bool=False
    else:
            for i in substrates.keys():
                    ub_list[i]=substrates[i]
            substrate_name=''
            substrate_value=0
            substrates_bool=True
    
    obj_name=product
    obj_target='maximize'
    biomass_id='BIOMASS_Ec_iML1515_core_75p37M'
    #pass_value={'A_tyr__L_e_con_1':tyr,'B_phpyr_e_con_2_reverse':php}
    EcoECM_FBA=Template_Concretemodel_3(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,reaction_list_C=reaction_list_C,reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,\
        obj_name=obj_name,obj_target=obj_target,set_obj_value=True,set_substrate_ini=True,substrate_name=substrate_name,\
        substrate_value=substrate_value,substrates=substrates,substrates_bool=substrates_bool,set_stoi_matrix=True,set_bound=True,\
        E_total=E_total,set_enzyme_constraint=False,set_part_enzyme_constraint=True,\
        biomass_id=biomass_id,set_A_biomass_ini=True,set_B_biomass_ini=True,set_C_biomass_ini=True,biomass_list={'A':bio_a,'B':bio_b,'C':bio_c},\
        set_pass_ini=False,set_pass_value=None,set_sub_ratio=False)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    #opt_ecm_pFBA.obj()
    return opt_ecm_FBA.obj()

