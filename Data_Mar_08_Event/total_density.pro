pro total_density

timespan,'2008-3-26/8',5,/h
thm_load_esa,probe='d',level='l2',datatype=['peef_density','peif_density']
thm_load_sst,probe='d',level='l2',datatype=['psef_density','psif_density']

get_data, 'thd_psef_density', data=d_psef
get_data, 'thd_psif_density', data=d_psif
get_data, 'thd_peef_density', data=d_peef
get_data, 'thd_peif_density', data=d_peif


d_psef_i={x:d_psif.x, y:interpol(d_psef.y,d_psef.x,d_psif.x)}
d_peef_i={x:d_psif.x, y:interpol(d_peef.y,d_peef.x,d_psif.x)}
d_peif_i={x:d_psif.x, y:interpol(d_peif.y,d_peif.x,d_psif.x)}

n_i={x:d_psif.x, y: d_psif.y + d_peif_i.y}
n_e={x:d_psif.x, y: d_psef_i.y + d_peef_i.y}

store_data,'thd_electron_density', data={x:n_e.x,y:n_e.y}
store_data,'thd_ion_density', data={x:n_i.x,y:n_i.y}

end
