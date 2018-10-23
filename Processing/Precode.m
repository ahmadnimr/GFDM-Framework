function xf =  Precode(D, params)
m_func = Modem_functions();
X = m_func.spread(D, params.K_set,params.EM1, params.EK2);
X = m_func.window(params.wtx,X, params.EW );
X = m_func.transform(X,'FD',params.EK3);
xf = m_func.allocation(X);
end
