function pin_legend_br(lg, ax)

    ax_pos = ax.Position;
    lg_w = lg.Position(3);
    lg.Position(1) = ax_pos(1) + ax_pos(3) - lg_w - 0.005;
    lg.Position(2) = ax_pos(2) + 0.005;

end
