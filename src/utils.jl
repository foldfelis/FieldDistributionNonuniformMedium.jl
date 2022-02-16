using Plots.PlotMeasures
using Plots

export
    plot_ϵ,
    plot_e_field

function plot_ϵ(s::Simulator; figsize=(600, 750), left_margin=-100px)
    ϵ = s.permittivity.ϵ

    return heatmap(
		axes(s.grid, 1), axes(s.grid, 2), ϵ',
		color=:algae,
		size=figsize, left_margin=left_margin, aspect_ratio=:equal
	)
end

function plot_e_field(s::Simulator; figsize=(600, 750), left_margin=-100px)
    ez = s.ez
    ϵ = s.permittivity.ϵ

    lim = maximum(abs.(ez))
    lim_ϵ = maximum(abs.(ϵ))

    p = plot(
        clim=(-lim, lim),  colorbar=false,
		size=figsize, left_margin=left_margin, aspect_ratio=:equal
    )
    p = heatmap!(
        p,
		axes(s.grid, 1), axes(s.grid, 2), ez',
		color=:coolwarm
	)
    p = contour!(
        p,
        axes(s.grid, 1), axes(s.grid, 2), lim .* ϵ' ./ lim_ϵ,
        color=:algae
    )

    return p
end
