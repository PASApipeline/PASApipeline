select distinct a.update_id, s.status_id, a.gene_id, a.model_id, a.alt_splice_flag, a.is_novel_flag from annotation_updates a, status s, status_link sl where a.is_valid = 1 and a.have_after = 1 and sl.compare_id = a.compare_id and s.status_id = sl.status_id and s.requires_update = 1 and sl.annot_update_id = a.update_id 

and s.status_id in (13,14,16) and a.compare_id = 3

order by s.status_id, a.update_id




