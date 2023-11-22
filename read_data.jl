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
    return stnid, name, state, lat, lon
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

function read_record(record_str::AbstractString)
    lines = split(record_str, "\n")

    # Parse metadata
    stnid, name, state, lat, lon = parse_metadata(lines[1])

    # Parse data lines
    observations = [parse_data(line) for line in lines[2:end]]

    return SingleDurationAnnMaxRecord(stnid, name, state, lon, lat, observations)
end

function read_all_records(filename::String)
    # Read the entire file into memory
    content = read(filename, String)

    # Split the content into individual records using a regular expression
    record_strs = split(content, r"\n{2,}")
    record_strs = [s for s in record_strs[2:end] if !isempty(s)]

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

    stn_df = DataFrame(; stnid=s.stnid, name=s.name, state=s.state, lon=s.lon, lat=s.lat)

    return obs_df, stn_df
end

function parse_file(filename::AbstractString)
    records = read_all_records(filename)
    dfs = [DataFrame(r) for r in records]
    obs_df = vcat([df[1] for df in dfs]...)
    stn_df = vcat([df[2] for df in dfs]...)
    return stn_df, obs_df
end
