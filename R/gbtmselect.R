#' gbtmselect
#'
#' As we all know, SAS "Traj" program can produce lots of files since we need to iterate our models. You can use this package to easily address the outcomes can converge them.
#' What you should do is nameing you outcome files as "oe/of/os_group_polynomial degree". Groups may be 1 to 5, and polynomial degree depends on the group numbers such as 33.
#' "oe" means parameter and covariance matrix estimates. "of" means group assignments and membership probabilities. "os" means parameter estimates used by TRAJPLOT macro.
#' More information can be found in "https://www.andrew.cmu.edu/user/bjones/index.htm".
#'
#' @param files_oe string vectors containing paths of GBTM outcome "oe" files.
#' @param files_of string vectors containing paths of GBTM outcome "of" files.
#' @param files_os string vectors containing paths of GBTM outcome "os" files.

#' @examples
#' gbtmselect(files_oe,files_of,files_os)
#' files_oe=c("oe//oe_1_1.sas7bdat","oe//oe_1_2.sas7bdat","oe//oe_1_3.sas7bdat","oe//oe_2_11.sas7bdat","oe//oe_2_22.sas7bdat","oe//oe_2_33.sas7bdat","oe//oe_3_011.sas7bdat","oe//oe_3_221.sas7bdat","oe//oe_3_333.sas7bdat","oe//oe_4_0111.sas7bdat","oe//oe_4_2221.sas7bdat","oe//oe_4_2331.sas7bdat")
#' files_of=c("of//of_1_1.sas7bdat","of//of_1_2.sas7bdat","of//of_1_3.sas7bdat","of//of_2_11.sas7bdat","of//of_2_22.sas7bdat","of//of_2_33.sas7bdat","of//of_3_011.sas7bdat","of//of_3_221.sas7bdat","of//of_3_333.sas7bdat","of//of_4_0111.sas7bdat","of//of_4_2221.sas7bdat","of//of_4_2331.sas7bdat")
#' files_os=c("os//os_1_1.sas7bdat","os//os_1_2.sas7bdat","os//os_1_3.sas7bdat","os//os_2_11.sas7bdat","os//os_2_22.sas7bdat","os//os_2_33.sas7bdat","os//os_3_011.sas7bdat","os//os_3_221.sas7bdat","os//os_3_333.sas7bdat","os//os_4_0111.sas7bdat","os//os_4_2221.sas7bdat","os//os_4_2331.sas7bdat")

#' @export
gbtmselect=function(files_oe,files_of,files_os){

  library(tidyverse)

  read_sas_gbtm=function(oe){
    # oe="oe/oe_1_1.sas7bdat"
    a=oe %>% str_split_fixed("_",n=3)
    b=a[3] %>% str_split_fixed(".s",n=2)
    oe_table=haven::read_sas(oe) %>% mutate(ngroups=as.numeric(a[2])) %>% mutate(order=(b[1]))
    return(oe_table)
  }

  oe=map_dfr(files_oe, read_sas_gbtm)
  of=map_dfr(files_of, read_sas_gbtm)
  os=map_dfr(files_os, read_sas_gbtm)

  e=oe %>% select("_LOGLIK_","_BIC1_","_BIC2_","ngroups","order") %>% unique()
  e=e[,c(4:5,1:3)]

  order=e[["order"]]
  order_space=function(order){
    n=length(order)
    order_sp=c()
    for (i in 1:n){
      if(str_length(order[i])==1){
        order_sp[i]=order[i]
      }else if(str_length(order[i])==2){
        order_sp[i]=paste(substr(order[i], 1, 1),substr(order[i],2,nchar(order[i])), sep =" ")
      }else if(str_length(order[i])==3){
        order_sp[i]=paste(substr(order[i], 1, 1),substr(order[i],2,2),substr(order[i],3,nchar(order[i])), sep =" ")
      }else if(str_length(order[i])==4){
        order_sp[i]=paste(substr(order[i], 1, 1),substr(order[i],2,2),substr(order[i],3,3),substr(order[i],4,nchar(order[i])), sep =" ")
      }else if(str_length(order[i])==5){
        order_sp[i]=paste(substr(order[i], 1, 1),substr(order[i],2,2),substr(order[i],3,3),substr(order[i],4,4),substr(order[i],5,nchar(order[i])), sep =" ")
      }
    }
    return(order_sp)
  }
  order_sp=tibble(order_space(order))

  e_sp=cbind(e,order_sp) %>% select(1,6,3,4,5,2)
  e_sp[[3]]=round(e_sp[[3]],3)
  e_sp[[4]]=round(e_sp[[4]],3)
  e_sp[[5]]=round(e_sp[[5]],3)

  s=os %>% select("PI","ngroups","order" ) %>% select(-ngroups)
  s=s %>% mutate(PI=round(PI,2)) %>% group_by(order) %>% summarise(mean_posterior_probabilities=str_c(PI[!is.na(PI)],collapse = "/")) %>%
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)

  of_1=of %>% group_by(order) %>% filter(GROUP==1) %>% select(order,GRP1PRB) %>% summarise(round(mean(GRP1PRB),2)) %>%
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)
  of_2=of %>% group_by(order) %>% filter(GROUP==2) %>% select(order,GRP2PRB) %>% summarise(round(mean(GRP2PRB),2)) %>%
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)
  of_3=of %>% group_by(order) %>% filter(GROUP==3) %>% select(order,GRP3PRB) %>% summarise(round(mean(GRP3PRB),2)) %>%
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)
  of_4=of %>% group_by(order) %>% filter(GROUP==4) %>% select(order,GRP4PRB) %>% summarise(round(mean(GRP4PRB),2)) %>%
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)
  of_5=of %>% group_by(order) %>% filter(GROUP==5) %>% select(order,GRP5PRB) %>% summarise(round(mean(GRP5PRB),2)) %>% 
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)

  of_sum=of_1 %>% left_join(of_2,by="order") %>% left_join(of_3,by="order") %>% left_join(of_4,by="order") %>% left_join(of_5,by="order")

  of_sum_1=of_sum %>% select(1,2) %>% mutate(per=`round(mean(GRP1PRB), 2)`) %>% select(-2)
  of_sum_2=of_sum %>% select(1,3) %>% mutate(per=`round(mean(GRP2PRB), 2)`) %>% select(-2)
  of_sum_3=of_sum %>% select(1,4) %>% mutate(per=`round(mean(GRP3PRB), 2)`) %>% select(-2)
  of_sum_4=of_sum %>% select(1,5) %>% mutate(per=`round(mean(GRP4PRB), 2)`) %>% select(-2)
  of_sum_5=of_sum %>% select(1,6) %>% mutate(per=`round(mean(GRP5PRB), 2)`) %>% select(-2)
  of_sum=rbind(of_sum_1,of_sum_2,of_sum_3,of_sum_4,of_sum_5)
  of_sum=of_sum %>% filter(!is.na(per))

  of_sum=of_sum %>% group_by(order) %>% summarise(Group_percentage=str_c(per[!is.na(per)],collapse = "/")) %>%
    mutate(i=nchar(order)) %>% arrange(i) %>% select(-i)

  sum_table=e_sp %>% left_join(s,by="order") %>% left_join(of_sum,by="order") %>% select(-6)
  colnames(sum_table)=c("Group","Polynomial degree","Log-Lik","BIC1","BIC2","Group percentage","Mean posterior probabilities")
  return(sum_table)
}

# devtools::document()

