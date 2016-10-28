##################################################################### UNIPROT ID
file_uniprot=open("uniprot_id.txt","r").readlines()
id_file_uniprot=[]
gene_file_uniprot=[]

for u in file_uniprot:
    tab=u.split("\t")
    genes=tab[-1].split(" ")
    for t in range(0,len(genes)):
        id_file_uniprot.append(tab[0])
        gene_file_uniprot.append(genes[t])

    if "\n" in gene_file_uniprot[-1]:
        gene_file_uniprot[-1] =gene_file_uniprot[-1][:-1]

print len(gene_file_uniprot)
print len(id_file_uniprot)
#if "LYN" in gene_file_uniprot:
 #   print "bieeen"
#######################################################################
id=[]
uniprot=[]
name=[]
goid=[]
goname=[]
fichier=["t-cell_activation.txt" ]# , "myeloid_dendritic_cell_activation.txt" ,  "follicular_dendritic_cell_activation.txt",  "immune_response.txt"]
for n in fichier:
    file=open(n, "r").readlines()
    for f in file:
###################################### take Uniprot AN and remove it from f
        uniprot.append(f[10:16])
        f=f[16:]
        while f[0] == " " or f[0] == "\t":
            f=f[1:]
###################################### take gene name
        if "\t" not in f:
            id.append(f[:f.index(" ")])
            f=f[f.index(" ")+1:]
        elif "\t" in f and len(f[:f.index("\t")])>15: ## if there's no tab in the middle
            id.append(f[:f.index(" ")])
            f=f[f.index(" ")+1:]
        elif "\t" in f and len(f[:f.index("\t")])<15:
            id.append(f[:f.index("\t")])
            f=f[f.index("\t")+1:]
######################################### take the Gene protein name
        name.append(f[:f.index("GO:")-1])
        while name[-1][-1] == " " or name[-1][-1] == "\t":
            name[-1]=name[-1][:-1]
        while name[-1][0] == " " or name[-1][0] == "\t":
            name[-1]=name[-1][1:]
        f=f[f.index("GO:"):]
######################################### take the GO id
        goid.append(f[:10])
        if "\t" in f:
            f=f[f.index("\t"):]
        if "\t" not in f:
            f=f[f.index(" "):]
        while f[0] == " " or f[0] == "\t":
            f=f[1:]
        goname.append(f[:-1])
#        if "response" in f[:-1] and "immune" not in f:
 #           print f[:-1]

############# CHECK BLANK SPACES
def blank(variable):
    for u in range(0, len(variable)):
        if len(variable[u]>2):
            while variable[u][0] == " ":
                variable[u]=variable[u][1:]
 ########################################################################

def VEGF(txt):
    cuenta_nogenes=0
    input=open(txt, "r").readlines()
    output=open(txt[:-4]+"_tcell.txt", "w")
    output.write("Metastatic-gene name\tMetastatic-Protein name\tMetastatic-GO Id\tMetastatic-GO name\tMetastatic-Value\tPrimary-gene name\tPrimary- Protein name\tPrimary-GO Id\tPrimary-GO name\tPrimary-Valu\
e\n")
    for f in input:
        #if "?" in f or "SLC35E2" in f:
            #break
        valuem= ""
        idm = "malm"
        genem= ""
        goidm= "NA"
        gonamem= ""
        absvaluem="NA"
        idp = "malp"
        genep= ""
        goidp= "NA"
        gonamep= ""
        valuep =""
        absvaluep="NA"
         if "|" in f and "tatic" not in f or "Gene name" not in f:
            f=f[:-1]
            tabs=f.split("\t")
            if len(tabs)==4:
                idm=tabs[0][:tabs[0].find("|")]
                if "NA" in tabs[1]:
                    valuem=0
                elif tabs[1] == "1":
                    valuem=1
            #                        absvaluem=1
                elif "-0." in tabs[1]:
                    valuem=tabs[1][tabs[1].find("0.")-1:]
                elif "." in tabs[1]:
                    valuem=tabs[1][tabs[1].find(".")-1:]
                #                        absvalem=abs(value)
#                if idm not in gene_file_uniprot:
 #                   cuenta_nogenes+=1
                if idm in gene_file_uniprot:
                    indicem=gene_file_uniprot.index(idm) ##index, numero
                    if id_file_uniprot[indicem] in uniprot:
                        indicem2=uniprot.index(id_file_uniprot[indicem])
                        genem=name[indicem2]
                        goidm=goid[indicem2]
                        gonamem=goname[indicem2]

                    elif gene_file_uniprot[indicem] in id:
                        indicem2=id.index(gene_file_uniprot[indicem])
                        genem=name[indicem2]
                        goidm=goid[indicem2]
                        gonamem=goname[indicem2]
################################################# PRIMARY
                if "|" in tabs[2]:
                    idp=tabs[2][:tabs[2].find("|")]
                if tabs[3] == "1":
                    valuep=1
 #                       absvaluep=1
                elif "-0." in tabs[3]:
                    valuep=tabs[3][tabs[3].find("0.")-1:]
                elif "." in tabs[3]:
                    valuep=tabs[3][tabs[3].find(".")-1:]
  #                      absvalep=abs(float(valuem))
                elif "NA" in tabs[3]:
                    valuep=0
                if idp in gene_file_uniprot:    #if our gene name is in the Uniprot list of gene names
                     indice=gene_file_uniprot.index(idp) ##index, numero
                     if id_file_uniprot[indice] in uniprot: #buscamos los AN correspondientes en uniprot que la lista de GO
                         indice2=uniprot.index(id_file_uniprot[indice])
                         genep=name[indice2]
                         goidp=goid[indice2]
                         gonamep=goname[indice2]
                     elif gene_file_uniprot[indice] in id:
                        indice2=id.index(gene_file_uniprot[indice])
                        genep=name[indice2]
                        goidp=goid[indice2]
                        gonamep=goname[indice2]

        output.write(idm+"\t"+genem+"\t"+goidm+"\t"+gonamem+"\t"+str(valuem)+"\t"+idp+"\t"+genep+"\t"+goidp+"\t"+gonamep+"\t"+str(valuep)+"\n")
    print cuenta_nogenes
print "C spearman"
VEGF("VEGFC_spearman.txt")
#print "C rsquared"
#VEGF("VEGFC_rsquared.txt")
print "A spearman"
VEGF("VEGFA_spearman.txt")
#print "A rsquared"
#VEGF("VEGFA_rsquared.txt")
print "D spearman"
VEGF("VEGFD_spearman.txt")
#print "D rsquared"
#VEGF("VEGFD_rsquared.txt")
