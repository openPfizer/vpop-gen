function dxdt = calc_dxdt_vandepas(t,x,p)

%% Species:
ifc = x(1); % [mmol]
ice = x(2); % [mmol]
hfc = x(3); % [mmol]
hce = x(4); % [mmol]
hdlfc = x(5); % [mmol]
hdlce = x(6); % [mmol]
nonhdl = x(7); % [mmol]
pfc = x(8); % [mmol]

%% Fluxes:
flux_hfc_synth = p(1); % Hepatic FC synthesis
flux_pfc_synth = p(2); % Peripheral FC synthesis
flux_ifc_synth = p(3); % Intestinal FC synthesis
flux_ifc_diet = p(4); % Intestinal FC from diet
flux_nonhdl_to_hfc = p(5)*nonhdl; % nonHDL (plasma) --> HFC
flux_hce_to_nonhdl = p(6)*hce; % HCE --> nonHDL (plasma)
flux_nonhdl_to_pfc = p(7)*nonhdl; % nonHDL (plasma) --> PFC
flux_pfc_to_hdlfc = p(8)*pfc; % RCT
flux_hdlfc_to_hdlce = p(9)*hdlfc; % Esterification on HDL
flux_hdlce_to_hdlfc = p(10)*hdlce; % Hydrolysis on HDL
flux_ice_to_nonhdl = p(11)*ice; % ICE --> nonHDL
flux_pfc_degrad = p(12)*pfc; % PFC --> *
flux_hdlfc_to_hfc = p(13)*hdlfc; % HDLfc --> HFC
flux_hfc_to_ifc = p(14)*hfc; % HFC --> IFC
flux_ifc_degrad = p(15)*ifc; % IFC --> *
flux_ifc_to_hdlfc = p(16)*ifc; % IFC --> HDLfc (plasma)
flux_hfc_to_hdlfc = p(17)*hfc; % HFC --> HDLfc (plasma)
flux_hfc_degrad = p(18)*hfc; % HFC --> *
flux_hfc_to_hce = p(19)*hfc; % HFC --> HCE
flux_ifc_to_ice = p(20)*ifc; % IFC --> ICE
flux_hdlce_to_nonhdlce = p(21)*hdlce; % HDLce --> nonHDLce

%% Species Balances:
dxdt_ifc = 
dxdt_ice = 
dxdt_hfc = 
dxdt_hce = 
dxdt_hdlfc = 
dxdt_hdlce = 
dxdt_nonhdl = 
dxdt_pfc = 

dxdt = [dxdt_ifc;dxdt_ice;dxdt_hfc;dxdt_hce;
    dxdt_hdlfc;dxdt_hdlce;dxdt_nonhdl;dxdt_pfc];

end