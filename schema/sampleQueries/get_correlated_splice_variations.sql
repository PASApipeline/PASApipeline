select distinct c.annotdb_asmbl_id, sv1.type, sv1.lend, sv1.rend, sv1.orient, sv2.type, sv2.lend, sv2.rend
from clusters c, cluster_link cl, splice_variation sv1, splice_variation sv2, alt_splice_link asl
where c.cluster_id = cl.cluster_id and cl.is_assembly = 1 and cl.cdna_acc = sv1.cdna_acc and sv1.sv_id = asl.sv_id_A and asl.sv_id_B = sv2.sv_id 
and sv1.sv_id < sv2.sv_id
