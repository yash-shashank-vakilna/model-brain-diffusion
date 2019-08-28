function [log_ind] = log_convert(lt)
b0 = 0.010000;
b1 = 0.010784;
log_ind = floor(log(lt/b0)/b1)+1;
end