import pandas as pd
import matplotlib.pyplot as plt
import math


def read_gauge(datapath, data: dict, limit=False, timeshift_true=False, trim_rows=0):
    """This parse function is capable of reading data produced by
    standard strain gauges by the tensile test machine in the
    composite lab."""
    df = pd.read_csv(datapath, encoding="latin1", delimiter="\t",
                     index_col=False)
    cross_area = data["cross_area"]
    if data["swap_strains"]:
        df.columns = ["time", "strainx", "strainy", "displacement", "load"]
    else:
        df.columns = ["time", "strainy", "strainx", "displacement", "load"]

    # Shifting the time, in case there something wrong with the data
    if timeshift_true:
        df["time"] = df["time"] + data["timeshift"]

    # Removing the part of the data gets a time less than zero
    df = df[df["time"] > 0]

    if limit is True:
        df = df[df["time"] < data["timelimit"]]
    elif limit > 0:
        df = df[df["time"] < limit]

    # for colname in df.columns:
    #     df[colname] = df[colname].astype(float)
    #     if df[colname].sum() < 0:
    #         # print(f"{colname=} had negative sum {df[colname].sum()}")
    #         df[colname] = df[colname] * -1

    df["strainx"] = df["strainx"] / 1000_000
    df["strainy"] = df["strainy"] / 1000_000
    df["stress"] = df["load"] / cross_area * 1000
    return df

def read_tensile3(datapath, metadata: dict, limit=False, length_adjust=0, timeshift_true=False, trim_rows=0):
    """This parse function is capable of the load - displacement data from the
    tensile test machine in the composite lab."""
    df = pd.read_csv(datapath, sep="\t").reset_index().drop([0, 1, 2])
    cross_area = metadata["width"]*metadata["thickness"]
    length = metadata["length"] + length_adjust

    # Name Columns appropriately
    df.columns = ["tull", "time", "displacement", "load"]
    df.drop(["tull"], 1)
    # df.columns = ["time", "load", "displacement"]
    # Convert values to floats
    for colname in df.columns[1:]:
        df[colname] = df[colname].str.replace(",", ".")
        df[colname] = df[colname].astype(float)

    if timeshift_true:
        df["time"] = df["time"] + metadata["timeshift"]

    if limit:
        df = df[df["time"] < metadata["timelimit"]]
    df["stress"] = df["load"] / cross_area
    df["strain"] = df["displacement"] / length

    # Dropping all rows that contain a NaN value
    df = df.dropna()

    return df

def read_tensile2(datapath, metadata: dict, limit=False, length_adjust=0, timeshift_true=False, trim_rows=0):
    """This parse function is capable of the load - displacement data from the
    tensile test machine in the composite lab."""
    df = pd.read_csv(datapath, sep=",").reset_index().drop([0, 1, 2])
    cross_area = metadata["width"]*metadata["thickness"]
    length = metadata["length"] + length_adjust

    # Name Columns appropriately
    df.columns = ["tull", "time", "load", "displacement", "also_extension"]
    df = df.drop(["tull", "also_extension"], 1)
    df = df.dropna()

    # df.columns = ["time", "load", "displacement"]
    # Convert values to floats
    #for colname in df.columns[1:]:
    #    df[colname] = df[colname].str.replace(",", ".")
    #    df[colname] = df[colname].astype(float)

    if timeshift_true:
        timeshift = metadata["timeshift"]
        df["time"] = df["time"] + metadata["timeshift"]


    # Then remove all lines that are negative time
    df = df[df["time"] >= 0]


    # Then adjust all loads and displacements with the first cell value
    starting_disp = df["displacement"].iloc[0]
    #print("starting displacement", starting_disp)
    #starting_load = df["load"].iloc[0]
    df["displacement"] = df["displacement"] - starting_disp

    if limit:
        df = df[df["time"] < metadata["timelimit"]]
    df["stress"] = df["load"] / cross_area
    df["strain"] = df["displacement"] / length

    # Dropping all rows that contain a NaN value
    df = df.dropna()

    return df

def read_tensile2_lbf(datapath, metadata: dict, limit=False, length_adjust=0, timeshift_true=False, trim_rows=0):
    """This parse function is capable of the load - displacement data from the
    tensile test machine in the composite lab."""
    df = pd.read_csv(datapath, sep=",").reset_index().drop([0, 1, 2])
    cross_area = metadata["width"]*metadata["thickness"]
    length = metadata["length"] + length_adjust

    # Name Columns appropriately
    df.columns = ["tull", "time", "load", "displacement", "also_extension"]
    df = df.drop(["tull", "also_extension"], 1)
    df = df.dropna()

    # adjust lbf to newtons
    df["load"] = df["load"]*4.448

    # df.columns = ["time", "load", "displacement"]
    # Convert values to floats
    #for colname in df.columns[1:]:
    #    df[colname] = df[colname].str.replace(",", ".")
    #    df[colname] = df[colname].astype(float)

    if timeshift_true:
        timeshift = metadata["timeshift"]
        df["time"] = df["time"] + metadata["timeshift"]


    # Then remove all lines that are negative time
    df = df[df["time"] >= 0]


    # Then adjust all loads and displacements with the first cell value
    starting_disp = df["displacement"].iloc[0]
    #print("starting displacement", starting_disp)
    #starting_load = df["load"].iloc[0]
    df["displacement"] = df["displacement"] - starting_disp

    if limit:
        df = df[df["time"] < metadata["timelimit"]]
    df["stress"] = df["load"] / cross_area
    df["strain"] = df["displacement"] / length

    # Dropping all rows that contain a NaN value
    df = df.dropna()

    return df

def read_flex_rectangular(datapath, metadata: dict, limit=False, length_adjust=0, trim_rows=0, timeshift_true=False):
    """This parse function is capable of the load - displacement data from the
    tensile test machine in the composite lab."""
    trim_rownums = list(range(0,trim_rows))
    df = pd.read_csv(datapath).reset_index().drop(trim_rownums)
    width = metadata["width"]
    thickness = metadata["thickness"]
    length = metadata["length"]
    cross_area = width*thickness
    length = metadata["length"] + length_adjust

    # Name Columns appropriately
    df.columns = ["what", "time", "load", "extra_displacement_column", "displacement"]

    # Remove the third column
    df = df[["time", "load", "displacement"]]

    # Convert values to floats
    for colname in df.columns:
        df[colname] = df[colname].astype(float)

    if timeshift_true:
        df["time"] = df["time"] + metadata["timeshift"]

    # Remove everything before 0
    df = df[df["time"] >= 0]

    if limit:
        df = df[df["time"] < metadata["timelimit"]]
        displacement_offset = df["displacement"].iloc[0]
        df["displacement"] = df["displacement"] - displacement_offset
    df["stress"] = 3*df["load"] *length / (2*width*thickness**2)
    df["strain"] = 6*df["displacement"]*thickness/(length**2)

    # Dropping all rows that contain a NaN value
    # df = df.dropna()
    return df

def read_flex_rectangular2(datapath, metadata: dict, limit=False, length_adjust=0, trim_rows=0, timeshift_true=False):
    """This parse function is capable of the load - displacement data from the
    tensile test machine in the composite lab."""
    trim_rownums = list(range(0,trim_rows))
    df = pd.read_csv(datapath).reset_index().drop(trim_rownums)
    width = metadata["width"]
    thickness = metadata["thickness"]
    length = metadata["length"]
    cross_area = width*thickness
    length = metadata["length"] + length_adjust

    # Name Columns appropriately
    df.columns = ["what", "time", "load", "displacement", "dontknow"]

    # Remove the third column
    df = df[["time", "load", "displacement"]]
    df = df.dropna()
    # Convert values to floats
    for colname in df.columns:
        df[colname] = df[colname].astype(float)

    if timeshift_true:
        df["time"] = df["time"] + metadata["timeshift"]

    # Remove everything before 0
    df = df[df["time"] >= 0]

    if limit:
        df = df[df["time"] < metadata["timelimit"]]
        displacement_offset = df["displacement"].iloc[0]
        df["displacement"] = df["displacement"] - displacement_offset
    df["stress"] = 3*df["load"] *length / (2*width*thickness**2)
    df["strain"] = 6*df["displacement"]*thickness/(length**2)

    # Dropping all rows that contain a NaN value
    # df = df.dropna()
    return df

def read_flex_circular(datapath, metadata: dict, limit=False, length_adjust=0, trim_rows=0, timeshift_true=False):
    """This parse function is capable of the load - displacement data from the
    tensile test machine in the composite lab."""
    trim_rownums = list(range(0,trim_rows))
    df = pd.read_csv(datapath).reset_index().drop(trim_rownums)
    radius = metadata["radius"]

    length = metadata["length"] + length_adjust

    # Name Columns appropriately
    df.columns = ["what", "time", "load", "extra_displacement_column", "displacement"]

    # Remove the third column
    df = df[["time", "load", "displacement"]]

    # Convert values to floats
    for colname in df.columns:
        df[colname] = df[colname].astype(float)

    if timeshift_true:
        df["time"] = df["time"] + metadata["timeshift"]

    # Remove everything before 0
    df = df[df["time"] >= 0]

    if limit:
        df = df[df["time"] < metadata["timelimit"]]
        displacement_offset = df["displacement"].iloc[0]
        df["displacement"] = df["displacement"] - displacement_offset
    df["stress"] = df["load"] * length / (math.pi * radius**3)
    df["strain"] = 6*df["displacement"]*2*radius/(length**2)

    # Dropping all rows that contain a NaN value
    df = df.dropna()

    return df


def chart_axes(xlabel="Strain", ylabel="Stress", title="Choose title"):
    fig = plt.figure()
    fig.title = title
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title)
    return axes

def get_modulus(strains, stresses, target_low_strain=0.001, target_high_strain=0.003):
    """This function returns an estimate for the elastic modulus based on the
    given strain locations. Example values for low and high strain is
    0.001 and 0.003."""
    strains = list(strains)
    stresses = list(stresses)
    low_stress = 0
    low_strain = 0
    high_stress = 0
    high_strain = 0
    for i, strain in enumerate(strains):
        if strain < target_low_strain:
            low_stress = stresses[i]
            low_strain = strain
        elif strain == target_low_strain:
            low_stress = stresses[i]
            low_strain = strain
            break
        else:
            break  # This means we already have the closest value lower than target

    for i, strain in enumerate(strains):
        if strain < target_high_strain:
            high_stress = stresses[i]
            high_strain = strain
        elif strain == target_high_strain:
            high_stress = stresses[i]
            high_strain = strain
            break
        else:
            break  # This means we already have the closest value lower than target

    # print(f"{low_stress=}, {high_stress=}, {low_strain=}, {high_strain=}")

    # Calculating modulus
    print("high strain", high_strain, "low strain", low_strain)
    modulus = (high_stress - low_stress) / (high_strain - low_strain)
    # print(f"Estimated high stress in target: {high_strain * modulus}")
    return modulus

def get_poissont(long_strains, trans_strains, target_low_strain, target_high_strain):
    """This function returns an estimate for the elastic modulus based on the
    given strain locations. For example 0.001 and 0.003."""
    long_strains = list(long_strains)
    trans_strains = list(trans_strains*-1)
    low_trans_strain = 0
    low_long_strain = 0
    high_trans_strain = 0
    high_long_strain = 0
    for i, strain in enumerate(long_strains):
        if strain < target_low_strain:
            low_trans_strain = trans_strains[i]
            low_long_strain = strain
        elif strain == target_low_strain:
            low_trans_strain = trans_strains[i]
            low_long_strain = strain
            break
        else:
            break  # This means we already have the closest value lower than target

    for i, strain in enumerate(long_strains):
        if strain < target_high_strain:
            high_trans_strain = trans_strains[i]
            high_long_strain = strain
        elif strain == target_high_strain:
            high_trans_strain = trans_strains[i]
            high_long_strain = strain
            break
        else:
            break  # This means we already have the closest value lower than target

    # print(f"{low_stress=}, {high_stress=}, {low_strain=}, {high_strain=}")

    # print(f"{high_trans_strain=}, {low_trans_strain=}, {high_long_strain=}, {low_long_strain=}")
    poissont_ratio = (high_trans_strain - low_trans_strain) / (high_long_strain - low_long_strain)
    # print(f"Estimated high stress in target: {high_strain * modulus}")
    return poissont_ratio


def get_strength():
    pass