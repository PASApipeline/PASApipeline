select distinct gene_id, model_id 
from annotation_link al, status_link sl, status s, cdna_link cdl, cluster_link cl 
where al.cdna_acc = sl.cdna_acc 
and sl.status_id = s.status_id 
and sl.compare_id = al.compare_id
and sl.cdna_acc = cdl.cdna_acc 
and cdl.cdna_acc = cl.cdna_acc 
and cl.is_assembly = 1 
and (s.requires_update = 1 or (s.requires_update = 0 and s.fails_incorporation = 0)) 
and cl.is_fli = 1
and al.compare_id = 3
