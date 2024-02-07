function [v_star, w_star, gv_star, gw_star, F_star] = ZKWolfe_try(v_0, w_0, dkv_0, dkw_0, pho, cstro, L, delta)
[F_0, gv_0, gw_0,] = ZK_CVCGFuncGrad_try(v_0, w_0);
gk_0 = transpose([gv_0(:); transpose(gw_0)]);
dk_0 = transpose([dkv_0(:); transpose(dkw_0)]);

r = 0;
gamma = (1 - cstro) / L;
eta_r = (pho^r) * gamma;
v_star = v_0 + eta_r * dkv_0; w_star = w_0 + eta_r * dkw_0;
[F_star, gv_star, gw_star] = ZK_CVCGFuncGrad_try(v_star, w_star);

t11 = real(gk_0 * dk_0');
t2_you = 0.9 * real(gk_0 * dk_0');

stop =0;
s = 0;
j = 0;


t1 = F_0 + delta * eta_r * t11;

gk_star = transpose([gv_star(:); transpose(gw_star)]);
t2_zuo = real(gk_star * dk_0');
if F_star <= t1 && t2_zuo >= t2_you
    return;
else
    while F_star > t1
        eta_r = (pho^j) * gamma;
        v_star = v_0 + eta_r * dkv_0; w_star = w_0 + eta_r * dkw_0;
        [F_star, gv_star, gw_star] = ZK_CVCGFuncGrad_try(v_star, w_star);
        t1 =  F_0 + delta * eta_r * t11;
        j = j + 1;
    end
    F_temp = F_star; gv_temp = gv_star; gw_temp = gw_star; v_temp = v_star; w_temp = w_star;
    while(1)
        gk_star = transpose([gv_star(:); transpose(gw_star)]);
        t2_zuo = real(gk_star * dk_0');
        if t2_zuo >= t2_you
            return
        else
            for s = 1:10
                eta_r = (eta_r + 2^s * (eta_r / 0.6 - eta_r)) * gamma;
                v_star = v_0 + eta_r * dkv_0; w_star = w_0 + eta_r * dkw_0;
                [F_star, gv_star, gw_star] = ZK_CVCGFuncGrad_try(v_star, w_star);
                gk_star = transpose([gv_star(:); transpose(gw_star)]);
                t2_zuo = real(gk_star * dk_0');
                if t2_zuo >= t2_you
                    t1 = F_0 + delta * eta_r * t11;
                    if F_star <= t1
                        return
                    end
                end
            end
        end
        break
    end
    F_star = F_temp; gv_star = gv_temp; gw_star = gw_temp; v_star = v_temp; w_star = w_temp;
end
end
