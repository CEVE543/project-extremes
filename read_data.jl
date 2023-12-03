using DataFrames
using Dates
using Unitful

struct Observation{T<:AbstractFloat}
    year::Int
    date::Date
    value::Quantity{T,Unitful.𝐋}
end

struct SingleDurationAnnMaxRecord{T<:AbstractFloat}
    stnid::AbstractString
    name::AbstractString
    state::AbstractString
    lon::T
    lat::T
    elevation::T
    records::Vector{Observation{T}}
end

year(s::SingleDurationAnnMaxRecord) = [obs.year for obs in s.records]
date(s::SingleDurationAnnMaxRecord) = [obs.date for obs in s.records]
value(s::SingleDurationAnnMaxRecord) = [obs.value for obs in s.records]

"""Parse the metadata from a given line"""
function parse_metadata(line::AbstractString)
    stnid = strip(line[1:7])
    name = strip(line[10:49])
    state = strip(line[52:53])
    lat = parse(Float64, strip(line[56:63]))

    # longitudes are a bit messier
    lon_chars = line[65:74]
    lon_chars = split(lon_chars, ",")[1]
    lon = parse(Float64, strip(lon_chars))
    elevation = parse(Float64, strip(line[77:end]))
    return stnid, name, state, lat, lon, elevation
end

"""Parse the data from a non-metadata line"""
function parse_data(line::AbstractString)
    date_str = strip(line[1:10])
    value_str = strip(line[12:end])
    date = Dates.Date(date_str, "mm/dd/yyyy")
    value = parse(Float64, value_str)u"inch"
    year = Dates.year(date)
    return Observation(year, date, value)
end

function read_record(lines::Vector{<:AbstractString})

    # Parse metadata
    stnid, name, state, lat, lon, elevation = parse_metadata(lines[1])

    # Parse data lines
    observations = [parse_data(line) for line in lines[2:end]]

    return SingleDurationAnnMaxRecord(stnid, name, state, lon, lat, elevation, observations)
end

function group_lines(lines)
    groups = Vector{String}[]  # List to store all groups
    current_group = String[]  # Current group of lines

    for line in lines
        if isempty(line)
            # When an empty line is encountered, store the current group (if it's not empty) and reset it
            if !isempty(current_group)
                push!(groups, current_group)
                current_group = String[]
            end
        else
            # Add non-empty lines to the current group
            push!(current_group, line)
        end
    end

    # Add the last group if it's not empty
    if !isempty(current_group)
        push!(groups, current_group)
    end

    return groups
end

function read_all_records(filename::String)
    # Read the entire file into memory
    all_lines = readlines(filename)[2:end]

    # Split the content into individual records using a regular expression
    record_strs = group_lines(all_lines)

    # Skip the first line and process each record, filtering out empty records
    records = [read_record(r) for r in record_strs]

    return records
end

function DataFrames.DataFrame(s::SingleDurationAnnMaxRecord)
    # Extract data for each observation
    years = [obs.year for obs in s.records]
    dates = [obs.date for obs in s.records]
    values = [obs.value for obs in s.records]

    # Create a DataFrame
    obs_df = DataFrame(;
        stnid=fill(s.stnid, length(s.records)), year=years, date=dates, value=values
    )

    stn_df = DataFrame(;
        stnid=s.stnid,
        name=s.name,
        state=s.state,
        lon=s.lon,
        lat=s.lat,
        elevation=s.elevation,
        nobs=length(s.records),
    )

    return obs_df, stn_df
end

function parse_file(filename::AbstractString)
    records = read_all_records(filename)
    dfs = [DataFrame(r) for r in records]
    obs_df = vcat([df[1] for df in dfs]...)
    stn_df = vcat([df[2] for df in dfs]...)
    return stn_df, obs_df
end
