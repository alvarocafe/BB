function calc_dfforma(qsi, eta)

dNdqsi = (1./4.)*[-(1. - eta);
                  (1. - eta);
                  (1. + eta);
                 -(1. + eta)];

dNdeta = (1./4.)*[-(1. - qsi);
                 -(1. + qsi);
                  (1. + qsi);
                  (1. - qsi)];
return dNdqsi, dNdeta
end

