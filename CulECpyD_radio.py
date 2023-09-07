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

def new_strip(str):#字符串分隔，酶信息修改工具
    while(1):
        if len(str)!=0 and str[0]=='('and str[-1]==')':
            str=str[1:-1]
        else:
            break
    return str

def sbml_to_excel(modelname,output):#模型SBML文件转为excel文件
    from cobra.io.dict import model_to_dict
    from cobra.io import read_sbml_model
    import pandas as pd
    model=read_sbml_model(modelname)
    model_dict=model_to_dict(model,sort=False)
    writer = pd.ExcelWriter(output)
    pd.DataFrame(model_dict['metabolites']).to_excel(writer,'metabolites',index=False)
    pd.DataFrame(model_dict['genes']).to_excel(writer,'genes',index=False)
    model_df=pd.DataFrame(model_dict['reactions'])
    model_df['reaction_eq'] = model_df.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=False))
    model_df['reaction_eq_name'] = model_df.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=True))
    model_df.to_excel(writer,'reactions',index=False)
    writer.save()

def isoenzyme_split(model):#同工酶反应切割
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

def Get_Models_info(model):#获取双菌模型基本信息
    reaction_list=[]
    reaction_list_A=[]
    reaction_list_B=[]
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
    return(reaction_list,reaction_list_A,reaction_list_B,metabolite_list,lb_list,ub_list,coef_matrix)

def Get_Concretemodel_Need_Data_json(model_file,reaction_kcat_MW_file,target):#以json格式文献为来源，将模型由json文件转为适用于pyomo处理的数学矩阵形式
    Concretemodel_Need_Data={}
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
    [reaction_list,reaction_list_A,reaction_list_B,metabolite_list,lb_list,ub_list,coef_matrix]=Get_Models_info(model)
    Concretemodel_Need_Data['reaction_list']=reaction_list
    Concretemodel_Need_Data['reaction_list_A']=reaction_list_A
    Concretemodel_Need_Data['reaction_list_B']=reaction_list_B
    Concretemodel_Need_Data['metabolite_list']=metabolite_list
    Concretemodel_Need_Data['lb_list']=lb_list
    Concretemodel_Need_Data['ub_list']=ub_list
    Concretemodel_Need_Data['coef_matrix']=coef_matrix
    return (Concretemodel_Need_Data)


def Model_Solve(model,solver):#标注运算器gurobi
    opt = pyo.SolverFactory(solver)
    opt.solve(model)
    return model

def unpdate_kcat_mw(reaction_list,kcat_mw,new_file):#更新模型对照酶动力学参数
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
    reaction_up_kact_MW.to_csv(new_file)
    return reaction_up_kact_MW

#Building

def exogenous_split(exogenous,meta_key):#切割途径外源引入反应信息
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

def add_new_ex(model,meta):#为原始模型新增外源反应
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

def get_new_model(model_name,exogenous,out_files,meta_key,add_product=False,product=None):#构建基于合成途径的子模型
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
    sbml_to_excel(out_files+'.xml',out_files+'.xlsx')  

def get_new_system(model_name,exogenous,model_A,model_B,meta_key,product):#联合子模型，构建联合系统
    (exogenous_A,exogenous_B)=exogenous_split(exogenous,meta_key)
    #A
    get_new_model(model_name,exogenous_A,model_A,meta_key)
    #B
    get_new_model(model_name,exogenous_B,model_B,meta_key,True,product)

def get_new_AB_model(model_origin_name,model_a_name,model_b_name,meta_test,DM_product,model_ab_json):#构建完整的基本双菌网络模型
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

def cut_compartment(id,compartment):#反应、代谢物id去除compartment后缀
    if id.endswith(compartment):
        id_cut=id[:-2]
    return id_cut

def EX_product(id):#生成目标产品外排反应
    return 'EX_'+cut_compartment(id,'_c')+'_e'

#simulation

def Template_Concretemodel_double(reaction_list=None,metabolite_list=None,coef_matrix=None,reaction_kcat_MW=None,lb_list=None,ub_list=None,reaction_list_A=None,reaction_list_B=None,\
    set_subs_ini=False,subs_list=None,subs_value=None,\
    set_biomass_ini=False,biomass_list=None,biomass_value=None,\
    set_bound=False,set_stoi_matrix=False,set_enzyme_constraint=False,set_part_enzyme_constraint=False,set_obj_single_E_value=False,E_total=None,\
    obj_name=None,obj_target=None,set_obj_value=False,set_target_ini=False,target_value=None,\
    set_obj_V_value=False,set_new_stoi_matrix=False,set_obj_E_value=False,\
    set_special_constraint=False,R_1=None,R_2=None):
    

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

    #To calculate the concentration variability of metabolites.
    # if set_obj_Met_value:   
    #     def set_obj_Met_value(m):
    #         return m.metabolite[obj_name]
    #     if obj_target=='maximize':
    #         Concretemodel.obj = Objective(rule=set_obj_Met_value, sense=maximize)
    #     elif obj_target=='minimize':
    #         Concretemodel.obj = Objective(rule=set_obj_Met_value, sense=minimize)

    #Adding flux balance constraints （FBA）
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
    #         return sum(coef_matrix[i,j]*m.reaction[j]  for j in reaction_list if (i,j) in coef_matrix.keys() )==0
    #     Concretemodel.set_new_stoi_matrix = Constraint( metabolite_list,rule=set_new_stoi_matrix)

    #Adding the upper and lower bound constraints of reaction flux
    if set_bound:
        def set_bound(m,j):
            return inequality(lb_list[j],m.reaction[j],ub_list[j])
        Concretemodel.set_bound = Constraint(reaction_list,rule=set_bound) 

    #Set the upper bound for substrate input reaction flux
    if set_subs_ini:
        def set_subs_ini(m,j): 
            return m.reaction[j] <= subs_value[j]
        Concretemodel.set_subs_ini = Constraint(subs_list,rule=set_subs_ini)   
    
    if set_biomass_ini:
        def set_biomass_ini(m,j): 
            return m.reaction[j] >= biomass_value[j]
        Concretemodel.set_biomass_ini = Constraint(biomass_list,rule=set_biomass_ini)  
    
    #Set the lower bound for target synthesis reaction flux
    if set_target_ini:
        def set_target_ini(m): 
            return m.reaction[obj_name] ==target_value
        Concretemodel.set_target_ini = Constraint(rule=set_target_ini)  

    #Adding enzymamic constraints
    if set_enzyme_constraint:
        def set_enzyme_constraint(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list)<= E_total
        Concretemodel.set_enzyme_constraint = Constraint(rule=set_enzyme_constraint)    
    
    #Adding Part A enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_A(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_A)<= E_total['A']
        Concretemodel.set_enzyme_constraint_A = Constraint(rule=set_enzyme_constraint_A)

    #Adding Part B enzymamic constraints
    if set_part_enzyme_constraint:
        def set_enzyme_constraint_B(m):
            return sum( m.reaction[j]/(reaction_kcat_MW.loc[j,'kcat_MW']) for j in reaction_kcat_MW.index if j in reaction_list_B)<= E_total['B']
        Concretemodel.set_enzyme_constraint_B = Constraint(rule=set_enzyme_constraint_B)

    if set_special_constraint:
        def set_special_constraint(m):
            return m.reaction[R_1] <= m.reaction[R_2]
        Concretemodel.set_special_constraint = Constraint(rule=set_special_constraint)


    return Concretemodel

def Model_Solve(model,solver):#标注运算器gurobi
    opt = pyo.SolverFactory(solver)
    opt.solve(model)
    return model

def ECM_FBA_subs(Concretemodel_Need_Data,subs,x,y,biomass,obj_name,E_total={'A':0.227,'B':0.227},set_special=False,R_1=None,R_2=None):#基于pyomo多底物的FBA分析方法
    # sub_id=[]
    # sub_value=[]
    # for i in subs.keys():
    #     sub_id.append(i)
    #     sub_value.append(subs[i])
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
    for i in subs.keys():
        #if subs[i] >= ub_list[i]:
        ub_list[i]=1000
    obj_name=obj_name
    obj_target='maximize'
    #biomass_id='BIOMASS_Ec_iML1515_core_75p37M'
    #pass_value={'A_tyr__L_e_con_1':mid_a,'B_phpyr_e_con_2_reverse':mid_b}
    EcoECM_FBA=Template_Concretemodel_double(reaction_list=reaction_list,metabolite_list=metabolite_list,coef_matrix=coef_matrix,\
        reaction_list_A=reaction_list_A,reaction_list_B=reaction_list_B,reaction_kcat_MW=reaction_kcat_MW,lb_list=lb_list,ub_list=ub_list,\
        obj_name=obj_name,obj_target=obj_target,set_obj_value=True,\
        set_subs_ini=True,subs_list=subs.keys(),subs_value=subs,\
        set_stoi_matrix=True,set_bound=True,E_total=E_total,\
        set_enzyme_constraint=False,set_part_enzyme_constraint=True,\
        set_biomass_ini=True,biomass_list=biomass.keys(),biomass_value=biomass,\
        set_special_constraint=set_special,R_1=R_1,R_2=R_2)
    opt_ecm_FBA=Model_Solve(EcoECM_FBA,'gurobi')
    opt_ecm_FBA.obj()
    return opt_ecm_FBA

def diff_subs_double_iML(Concretemodel_Need_Data,double,subs,obj_name):#
    (x,y)=double
    #bio value from ECMpy https://doi.org/10.3390/biom12010065 
    bio_a=0.45
    bio_b=0.65
    biomass_list={'A_BIOMASS_Ec_iML1515_core_75p37M':bio_a,'B_BIOMASS_Ec_iML1515_core_75p37M':bio_b}
    a=ECM_FBA_subs(Concretemodel_Need_Data,subs,x,y,biomass_list,obj_name,E_total={'A':0.227,'B':0.227},set_special=False,R_1=None,R_2=None)
    print('A:B='+str(x)+':'+str(y)+'→obj:')
    print(a.obj())
    return a.obj()

def diff_subs_iML_iCW(Concretemodel_Need_Data,double,subs,obj_name,set_special=False):#
    (x,y)=double
    #Setting ensures the growth state of the strain
    bio_a=0.1
    bio_b=0.1
    biomass_list={'A_BIOMASS_Ec_iML1515_core_75p37M':bio_a,'B_CG_biomass_cgl_ATCC13032':bio_b}
    a=ECM_FBA_subs(Concretemodel_Need_Data,subs,x,y,biomass_list,obj_name,{'A':0.227,'B':0.129},set_special=set_special,R_1='EX_15dap_e',R_2='EX_lys__L_e')
    print('A:B='+str(x)+':'+str(y)+'→obj:')
    print(a.obj())
    return a.obj()

def unpdate_kcat_mw_diff_hosts(reaction_list,kcat_file_A,kcat_file_B,new_file):
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

