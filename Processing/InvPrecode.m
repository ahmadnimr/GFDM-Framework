function D =  InvPrecode(y, params)
m_func = Modem_functions();
[K,M] = size(params.wrx);
D = m_func.deallocation(y,K,M,'FD');
D = m_func.detransform(D,'FD', params.EK3);
D = m_func.window(params.wrx, D, params.EW );
D = m_func.despread(D, params.K_set,params.EM1, params.EK2);
end