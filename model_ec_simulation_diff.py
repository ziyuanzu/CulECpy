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
import pyomo.environ as pyo
from pyomo.environ import *
from pyomo.opt import SolverFactory

def Template_Concretemodel_double_diff(reaction_list=None,metabolite_list=None,coef_matrix=None,reaction_kcat_MW=None,lb_list=None,ub_list=None,reaction_list_A=None,reaction_list_B=None,\
    set_substrate_ini=False,substrates=None,\
    set_A_biomass_ini=None,set_B_biomass_ini=None,biomass_value_A=None,biomass_value_B=None,biomass_id_A=None,biomass_id_B=None,\
    set_bound=False,set_stoi_matrix=False,\
    set_part_enzyme_constraint=False,set_obj_single_E_value=False,E_total_A=None,E_total_B=None,\
    set_obj_value=False,obj_name=None,obj_target=None,\
    set_target_ini=False,target_name=None,target_value=None,\
    set_obj_V_value=False,set_obj_E_value=False,x=None,y=None):
    

    Concretemodel = ConcreteModel()
    Concretemodel.reaction = pyo.Var(reaction_list,  within=NonNegativeReals)
    Concretemodel.z = pyo.Var(reaction_list,  within=pyo.Binary)
    Concretemodel.x = pyo.Var()
    Concretemodel.y = pyo.Var()

    #Set upper and lower bounds of metabolite concentration
    # if set_metabolite:
    #     def set_metabolite(m,i):
    #         return  inequality(metabolites_lnC.loc[i,'lnClb'], m.metabolite[i], metabolites_lnC.loc[i,'lnCub'])
    #     Concretemodel.set_metabolite= Constraint(metabolite_list,rule=set_metabolite)     

    #Minimizing the flux sum of pathway (pFBA)
    if set_obj_V_value:             
        def set_obj_V_value(m):
            return sum(m.reaction[j] for j in reaction_list)
        Concretemodel.obj = Objective(rule=set_obj_V_value, sense=minimize)  
    
    #Set the maximum flux as the object function
    if set_obj_value:   
        def set_obj_value(m):
            return m.reaction[obj_name]#+m.reaction['A_'+biomass_id_A]+m.reaction['B_'+biomass_id_B]
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

    #Adding flux balance constraints （FBA）
    if set_stoi_matrix:
        def set_stoi_matrix(m,i):
            return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
        Concretemodel.set_stoi_matrix = Constraint( metabolite_list,rule=set_stoi_matrix)

    #Adding the upper and lower bound constraints of reaction flux
    if set_bound:
        def set_bound(m,j):
            return inequality(lb_list[j],m.reaction[j],ub_list[j])
        Concretemodel.set_bound = Constraint(reaction_list,rule=set_bound) 
    
     #Set the upper bound for substrate input reaction flux
    if set_substrate_ini :
        def set_substrates_ini(m,j): 
            return m.reaction[j]<= substrates[j]
        Concretemodel.set_substrates_ini = Constraint(substrates.keys(),rule=set_substrates_ini)  
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_A_biomass_ini:
        def set_A_biomass_ini(m): 
            return m.reaction['A_'+biomass_id_A] >=biomass_value_A
        Concretemodel.set_A_biomass_ini = Constraint(rule=set_A_biomass_ini)  
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_B_biomass_ini:
        def set_B_biomass_ini(m): 
            return m.reaction['B_'+biomass_id_B] >=biomass_value_B
        Concretemodel.set_B_biomass_ini = Constraint(rule=set_B_biomass_ini)  
    
    #Set the lower bound for target synthesis reaction flux
    if set_target_ini:
        def set_target_ini(m): 
            return m.reaction[target_name] >=target_value
        Concretemodel.set_target_ini = Constraint(rule=set_target_ini)  
    
    #Adding Part A enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_A(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_A)<= E_total_A
        Concretemodel.set_enzyme_constraint_A = Constraint(rule=set_enzyme_constraint_A)

    #Adding Part B enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_B(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_B)<= E_total_B
        Concretemodel.set_enzyme_constraint_B = Constraint(rule=set_enzyme_constraint_B)
    

        # def special_constraint2(m):
        #     return m.reaction['A_'+biomass_id_A]*x+m.reaction['B_'+biomass_id_B]*y>=0.46
        # Concretemodel.special_constraint2 = Constraint(rule=special_constraint2)
    return Concretemodel

def Model_Solve(model,solver):#标注运算器gurobi
    opt = pyo.SolverFactory(solver)
    opt.solve(model)
    return model

def ECM_FBA_diff(Concretemodel_Need_Data,substrates,x,y,bio_a,bio_b,obj_name,obj_target='maximize',set_target_ini=False,target_name=None,target_value=None,special_constraint=False):#基于pyomo的FBA分析方法（glc）
    E_total_A=0.227
    E_total_B=0.227
    sub_a='EX_glc__D_e_reverse'
    #ratio=4
    reaction_list_A=deepcopy(Concretemodel_Need_Data['reaction_list_A'])
    reaction_list_B=deepcopy(Concretemodel_Need_Data['reaction_list_B'])
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
                        # if 'akg_e_con_' in j:
                #                 print ('B:'+str(coef_matrix[i,j]))
                #     print ('ij:'+str(coef_matrix[i,j]))
    reaction_kcat_MW=deepcopy(Concretemodel_Need_Data['reaction_kcat_MW'])
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    biomass_id_a='BIOMASS_Ec_iML1515_core_75p37M'
    biomass_id_b='BIOMASS_Ec_iML1515_core_75p37M'    
    #pass_value={'A_tyr__L_e_con_1':tyr,'B_phpyr_e_con_2_reverse':php}
    EcoECM_FBA=Template_Concretemodel_double_diff(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,\
        set_substrate_ini=True,substrates=substrates,\
        set_A_biomass_ini=True,set_B_biomass_ini=True,biomass_value_A=bio_a,biomass_value_B=bio_b,biomass_id_A=biomass_id_a,biomass_id_B=biomass_id_b,\
        set_bound=True,set_stoi_matrix=True,set_obj_single_E_value=False,\
        set_part_enzyme_constraint=True,E_total_A=E_total_A,E_total_B=E_total_B,\
        set_obj_value=True,obj_name=obj_name,obj_target=obj_target,\
        set_target_ini=set_target_ini,target_name=target_name,target_value=target_value,\
        set_obj_V_value=False,set_obj_E_value=False,x=x,y=y)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    #opt_ecm_pFBA.obj()
    return opt_ecm_FBA

def ECM_pFBA_diff(Concretemodel_Need_Data,substrates,x,y,bio_a,bio_b,obj_name,obj_target='minimize',set_target_ini=False,target_name=None,target_value=None,special_constraint=False):#基于pyomo的FBA分析方法（glc）
    E_total_A=0.227
    E_total_B=0.227
    sub_a='EX_glc__D_e_reverse'
    #ratio=4
    reaction_list_A=deepcopy(Concretemodel_Need_Data['reaction_list_A'])
    reaction_list_B=deepcopy(Concretemodel_Need_Data['reaction_list_B'])
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
                        # if 'akg_e_con_' in j:
                #                 print ('B:'+str(coef_matrix[i,j]))
                #     print ('ij:'+str(coef_matrix[i,j]))
    reaction_kcat_MW=deepcopy(Concretemodel_Need_Data['reaction_kcat_MW'])
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    biomass_id_a='BIOMASS_Ec_iML1515_core_75p37M'
    biomass_id_b='BIOMASS_Ec_iML1515_core_75p37M'    
    #pass_value={'A_tyr__L_e_con_1':tyr,'B_phpyr_e_con_2_reverse':php}
    EcoECM_FBA=Template_Concretemodel_double_diff(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,\
        set_substrate_ini=True,substrates=substrates,\
        set_A_biomass_ini=True,set_B_biomass_ini=True,biomass_value_A=bio_a,biomass_value_B=bio_b,biomass_id_A=biomass_id_a,biomass_id_B=biomass_id_b,\
        set_bound=True,set_stoi_matrix=True,set_obj_single_E_value=False,\
        set_part_enzyme_constraint=True,E_total_A=E_total_A,E_total_B=E_total_B,\
        set_obj_value=False,obj_name=obj_name,obj_target=obj_target,\
        set_target_ini=set_target_ini,target_name=target_name,target_value=target_value,\
        set_obj_V_value=True,set_obj_E_value=False,x=x,y=y)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    #opt_ecm_pFBA.obj()
    return opt_ecm_FBA

def FBA_diff(Concretemodel_Need_Data,substrates,x,y,bio_a,bio_b,obj_name,obj_target='maximize',set_target_ini=False,target_name=None,target_value=None,special_constraint=False):#基于pyomo的FBA分析方法（glc）
    E_total_A=0.227
    E_total_B=0.227
    sub_a='EX_glc__D_e_reverse'
    #ratio=4
    reaction_list_A=deepcopy(Concretemodel_Need_Data['reaction_list_A'])
    reaction_list_B=deepcopy(Concretemodel_Need_Data['reaction_list_B'])
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
                        # if 'akg_e_con_' in j:
                #                 print ('B:'+str(coef_matrix[i,j]))
                #     print ('ij:'+str(coef_matrix[i,j]))
    reaction_kcat_MW=deepcopy(Concretemodel_Need_Data['reaction_kcat_MW'])
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    biomass_id_a='BIOMASS_Ec_iML1515_core_75p37M'
    biomass_id_b='BIOMASS_Ec_iML1515_core_75p37M'    
    #pass_value={'A_tyr__L_e_con_1':tyr,'B_phpyr_e_con_2_reverse':php}
    EcoECM_FBA=Template_Concretemodel_double_diff(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,\
        set_substrate_ini=True,substrates=substrates,\
        set_A_biomass_ini=True,set_B_biomass_ini=True,biomass_value_A=bio_a,biomass_value_B=bio_b,biomass_id_A=biomass_id_a,biomass_id_B=biomass_id_b,\
        set_bound=True,set_stoi_matrix=True,set_obj_single_E_value=False,\
        set_part_enzyme_constraint=False,E_total_A=E_total_A,E_total_B=E_total_B,\
        set_obj_value=True,obj_name=obj_name,obj_target=obj_target,\
        set_target_ini=set_target_ini,target_name=target_name,target_value=target_value,\
        set_obj_V_value=False,set_obj_E_value=False,x=x,y=y)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    #opt_ecm_pFBA.obj()
    return opt_ecm_FBA

def pFBA_diff(Concretemodel_Need_Data,substrates,x,y,bio_a,bio_b,obj_name,obj_target='minimize',set_target_ini=False,target_name=None,target_value=None,special_constraint=False):#基于pyomo的FBA分析方法（glc）
    E_total_A=0.227
    E_total_B=0.227
    sub_a='EX_glc__D_e_reverse'
    #ratio=4
    reaction_list_A=deepcopy(Concretemodel_Need_Data['reaction_list_A'])
    reaction_list_B=deepcopy(Concretemodel_Need_Data['reaction_list_B'])
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
                        # if 'akg_e_con_' in j:
                #                 print ('B:'+str(coef_matrix[i,j]))
                #     print ('ij:'+str(coef_matrix[i,j]))
    reaction_kcat_MW=deepcopy(Concretemodel_Need_Data['reaction_kcat_MW'])
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    biomass_id_a='BIOMASS_Ec_iML1515_core_75p37M'
    biomass_id_b='BIOMASS_Ec_iML1515_core_75p37M'    
    #pass_value={'A_tyr__L_e_con_1':tyr,'B_phpyr_e_con_2_reverse':php}
    EcoECM_FBA=Template_Concretemodel_double_diff(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,\
        set_substrate_ini=True,substrates=substrates,\
        set_A_biomass_ini=True,set_B_biomass_ini=True,biomass_value_A=bio_a,biomass_value_B=bio_b,biomass_id_A=biomass_id_a,biomass_id_B=biomass_id_b,\
        set_bound=True,set_stoi_matrix=True,set_obj_single_E_value=False,\
        set_part_enzyme_constraint=False,E_total_A=E_total_A,E_total_B=E_total_B,\
        set_obj_value=False,obj_name=obj_name,obj_target=obj_target,\
        set_target_ini=set_target_ini,target_name=target_name,target_value=target_value,\
        set_obj_V_value=True,set_obj_E_value=False,x=x,y=y)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    #opt_ecm_pFBA.obj()
    return opt_ecm_FBA

def unpdate_kcat_mw_diff(reaction_list,kcat_file_A,kcat_file_B,new_file):
    reaction_up_kact_MW=pd.DataFrame()
    kcat_A=pd.read_csv(kcat_file_A,index_col=0)
    kcat_B=pd.read_csv(kcat_file_B,index_col=0)
    for i in reaction_list:
        if i.startswith('A_'):
            id=i[2:]
            if id in kcat_A.index:
                reaction_up_kact_MW.at[i,'kcat_MW']=kcat_A.at[id,'kcat_MW']
        if i.startswith('B_'):
            id=i[2:]
            if id in kcat_B.index:
                reaction_up_kact_MW.at[i,'kcat_MW']=kcat_B.at[id,'kcat_MW'] 
    reaction_up_kact_MW.to_csv(new_file)
    return reaction_up_kact_MW


def Template_Concretemodel_double_imlicw(reaction_list=None,metabolite_list=None,coef_matrix=None,reaction_kcat_MW=None,lb_list=None,ub_list=None,reaction_list_A=None,reaction_list_B=None,\
    set_substrate_ini=False,substrates=None,\
    set_A_biomass_ini=None,set_B_biomass_ini=None,biomass_value_A=None,biomass_value_B=None,biomass_id_A=None,biomass_id_B=None,\
    set_bound=False,set_stoi_matrix=False,\
    set_part_enzyme_constraint=False,set_obj_single_E_value=False,E_total_A=None,E_total_B=None,\
    set_obj_value=False,obj_name=None,obj_target=None,\
    set_target_ini=False,target_name=None,target_value=None,\
    set_obj_V_value=False,set_obj_E_value=False,special_constraint=False,x=None,y=None):
    

    Concretemodel = ConcreteModel()
    Concretemodel.reaction = pyo.Var(reaction_list,  within=NonNegativeReals)
    Concretemodel.z = pyo.Var(reaction_list,  within=pyo.Binary)
    Concretemodel.x = pyo.Var()
    Concretemodel.y = pyo.Var()

    #Set upper and lower bounds of metabolite concentration
    # if set_metabolite:
    #     def set_metabolite(m,i):
    #         return  inequality(metabolites_lnC.loc[i,'lnClb'], m.metabolite[i], metabolites_lnC.loc[i,'lnCub'])
    #     Concretemodel.set_metabolite= Constraint(metabolite_list,rule=set_metabolite)     

    #Minimizing the flux sum of pathway (pFBA)
    if set_obj_V_value:             
        def set_obj_V_value(m):
            return sum(m.reaction[j] for j in reaction_list)
        Concretemodel.obj = Objective(rule=set_obj_V_value, sense=minimize)  
    
    #Set the maximum flux as the object function
    if set_obj_value:   
        def set_obj_value(m):
            return m.reaction[obj_name]#+m.reaction['A_'+biomass_id_A]+m.reaction['B_'+biomass_id_B]
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

    #Adding flux balance constraints （FBA）
    if set_stoi_matrix:
        def set_stoi_matrix(m,i):
            return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
        Concretemodel.set_stoi_matrix = Constraint( metabolite_list,rule=set_stoi_matrix)

    #Adding the upper and lower bound constraints of reaction flux
    if set_bound:
        def set_bound(m,j):
            return inequality(lb_list[j],m.reaction[j],ub_list[j])
        Concretemodel.set_bound = Constraint(reaction_list,rule=set_bound) 
    
     #Set the upper bound for substrate input reaction flux
    if set_substrate_ini :
        def set_substrates_ini(m,j): 
            return m.reaction[j]<= substrates[j]
        Concretemodel.set_substrates_ini = Constraint(substrates.keys(),rule=set_substrates_ini)  
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_A_biomass_ini:
        def set_A_biomass_ini(m): 
            return m.reaction['A_'+biomass_id_A] >=biomass_value_A
        Concretemodel.set_A_biomass_ini = Constraint(rule=set_A_biomass_ini)  
        
    #Set the lower bound for biomass synthesis reaction flux
    if set_B_biomass_ini:
        def set_B_biomass_ini(m): 
            return m.reaction['B_'+biomass_id_B] >=biomass_value_B
        Concretemodel.set_B_biomass_ini = Constraint(rule=set_B_biomass_ini)  
    
    #Set the lower bound for target synthesis reaction flux
    if set_target_ini:
        def set_target_ini(m): 
            return m.reaction[target_name] >=target_value
        Concretemodel.set_target_ini = Constraint(rule=set_target_ini)  
    
    #Adding Part A enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_A(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_A)<= E_total_A
        Concretemodel.set_enzyme_constraint_A = Constraint(rule=set_enzyme_constraint_A)

    #Adding Part B enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_B(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_B)<= E_total_B
        Concretemodel.set_enzyme_constraint_B = Constraint(rule=set_enzyme_constraint_B)
    
    if special_constraint:
        def special_constraint(m):
            return m.reaction['A_pdx']==m.reaction['A_LYSDC_num1']+m.reaction['A_LYSDC_num2']
        Concretemodel.special_constraint = Constraint(rule=special_constraint)
        # def special_constraint2(m):
        #     return m.reaction['A_'+biomass_id_A]*x+m.reaction['B_'+biomass_id_B]*y>=0.46
        # Concretemodel.special_constraint2 = Constraint(rule=special_constraint2)
    return Concretemodel

def ECM_FBA_imlicw(Concretemodel_Need_Data,substrates,x,y,bio_a,bio_b,obj_name,obj_target='maximize',set_target_ini=False,target_name=None,target_value=None,special_constraint=False):#基于pyomo的FBA分析方法（glc）
    E_total_A=0.227
    E_total_B=0.129
    sub_a='EX_glc__D_e_reverse'
    sub_b='EX_xyl__D_e_reverse'
    ratio=4
    reaction_list_A=deepcopy(Concretemodel_Need_Data['reaction_list_A'])
    reaction_list_B=deepcopy(Concretemodel_Need_Data['reaction_list_B'])
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
                        # if 'akg_e_con_' in j:
                #                 print ('B:'+str(coef_matrix[i,j]))
                #     print ('ij:'+str(coef_matrix[i,j]))
    reaction_kcat_MW=deepcopy(Concretemodel_Need_Data['reaction_kcat_MW'])
    lb_list=Concretemodel_Need_Data['lb_list']
    ub_list=Concretemodel_Need_Data['ub_list']
    biomass_id_a='BIOMASS_Ec_iML1515_core_75p37M'
    biomass_id_b='CG_biomass_cgl_ATCC13032'    
    #pass_value={'A_tyr__L_e_con_1':tyr,'B_phpyr_e_con_2_reverse':php}
    EcoECM_FBA=Template_Concretemodel_double_imlicw(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,\
        set_substrate_ini=True,substrates=substrates,\
        set_A_biomass_ini=True,set_B_biomass_ini=True,biomass_value_A=bio_a,biomass_value_B=bio_b,biomass_id_A=biomass_id_a,biomass_id_B=biomass_id_b,\
        set_bound=True,set_stoi_matrix=True,set_obj_single_E_value=False,\
        set_part_enzyme_constraint=True,E_total_A=E_total_A,E_total_B=E_total_B,\
        set_obj_value=True,obj_name=obj_name,obj_target='maximize',\
        set_target_ini=set_target_ini,target_name=target_name,target_value=target_value,\
        set_obj_V_value=False,set_obj_E_value=False,special_constraint=special_constraint,x=x,y=y)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    #opt_ecm_pFBA.obj()
    return opt_ecm_FBA
