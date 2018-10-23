function D = Equalization(yf, eq, params)
% FD equalization
% yf: N x 1 received data
% Params: precoding parmeters
if params.EW % not OFDM
    u_func = Utility_functions();
    Vyf = u_func.V(yf, params.K, params.M);
    Vyf= eq.H_eq.*Vyf;
    yf_eq = u_func.iV(Vyf);
else % OFDM No Channel equlization is required
    yf_eq = yf;
end
D = InvPrecode(yf_eq, params); % inverse precoding.. See bellow
% Scale with the gain prior to demapping
D = D./eq.Ag;
end

