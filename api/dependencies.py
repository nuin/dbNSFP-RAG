"""Resolve API database providers at request time."""


def get_database(genome_build="GRCh37"):
    from api import server

    return server.get_database(genome_build)
