function repin_all(legs, axs)

    drawnow;

    for j = 1:length(axs)
        pin_legend_br(legs(j), axs(j));
    end

end
