function load_description_file(datadir, subdir, qname)
    qpath = joinpath(datadir, subdir, qname)
    df = CSV.read(qpath, DataFrame, header=1)

    # mass[kg] => g/mol * mass_factor
    mass_factor = 1.0e-3/6.02214076e23

    isotope_id = df[!,:localID]
    isotope_a  = df[!,:Abundance]
    isotope_m  = df[!,"MolarMass/g·mol-1"] .* mass_factor
    qpaths     = df[!,"Q(fullrange)"]
    gis        = df[!,:gi]

    isotope_id, isotope_a, isotope_m, gis, qpaths
end
function get_melcule_and_isotope_ids()
    @infoe @sprintf("max_isotope_id = %d", max_isotope_id)
    iso_ids = Vector{Int64}(undef,0)
    for iline in 1:nb_lines
        iso_id = spectral_data[12, iline]
        isotope_id = floor(Int64, iso_id)
        if isotope_id <= max_isotope_id
            push!(iso_ids, isotope_id)
        end
    end
end
