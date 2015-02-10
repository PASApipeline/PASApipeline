select c.annotdb_asmbl_id, cdl.cdna_acc, cdl.alignment from clusters c, cluster_link cl, cdna_link cdl where c.cluster_id = cl.cluster_id and cl.is_assembly = 1 and cl.cdna_acc = cdl.cdna_acc

