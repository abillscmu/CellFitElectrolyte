using Turing

function neg()
    return truncated(Normal(), upper=0)
end

function pos()
    return truncated(Normal(), lower=0)
end

@model function fit_linear(data_dict, cycle_1, cycle_2, cell, deg_df)
    #SEI NEGATIVE
    sei_neg_pore ~ neg()
    sei_neg_resist ~ pos()
    sei_neg_n_li ~ neg()
    #SEI POSITIVE
    sei_pos_pore ~ neg()
    sei_pos_resist ~ pos()
    sei_pos_n_li ~ neg()
    #PLATING NEGATIVE
    plating_pore ~ neg()
    plating_n_li ~ neg()
    #CLOGGING NEGATIVE


end

    


