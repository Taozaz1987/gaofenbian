function V_ALFMT_r=calculate_Vr(V_ALFMT_s,ALFMT_s_hs,phi)
miu=min(ALFMT_s_hs(:));
A_max=max(ALFMT_s_hs(:));
V_ALFMT_r=V_ALFMT_s./(ALFMT_s_hs+miu*A_max).*exp(-1i.*phi);
