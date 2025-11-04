function data_rmse = rmse(data, u)
    data_rmse = sqrt(mean((data-u).^2));
end