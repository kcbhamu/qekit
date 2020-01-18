module qekit
    using Shell
    using HDF5
    using Supressor
    using DelimitedFiles

    struct qe
        input  :: String
        output :: String
        prefix :: String
        fermi  :: Float64
        energy :: Float64
        magnet :: Dict
        charge :: Dict

        qe(input,output) = new(input,ouput,"",0,0,Dict(),Dict())
    end

    struct wannier90
    end

    getMagnet(crystal::qe) = crystal.magnet
    getFermi(crystal::qe) = crystal.fermi
    getEnergy(crystal::qe) = crystal.energy
    getCharge(crysta::qe) = crystal.charge

    function pw_band(crystal::qe)
    end

    function pw_dos(crystal::qe)
    end

    function pw_pdos(crystal::qe)
    end

    function get_w90Input(crystal::qe)
    end

    function w90_band(w90::wannier90)
    end

    function w90_pdos(w90::wannier90)
    end

    function w90_boltz(w90::wannier90)
    end

    function w90_berry(w90::wannier90)
    end

    function qekit_info()
    end

    function w90_lib()
    end

    function pw_lib()
    end

end
