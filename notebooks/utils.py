#!/usr/bin/env python
# coding: utf-8


def join_string_columns(df, column, join_cols, keep=False, delimiter="|"):
    import pandas as pd

    """
    Join multiple string columns.

    This function joins multiple text columns using the specified delimiter.
    Empty strings '' and None will be ignored.

    Parameters
    ----------
    df : DataFrame
        First number to add.
    column : str
        Name of column with the joined string.
    join_cols : list
        List of column names.
    delete : bool, default False
        If True, delete join_cols, but exclude result column if within the join_cols list.
    delimiter : str, default '|'
        Delimiter used to join columns

    Returns
    -------
    DataFrame
        Dataframe with column of the concatenated string.


    Examples
    --------
    >>> df = pd.DataFrame([{'A': 'a', 'B': 'b', 'C': 'c'},
    >>>                    {'A': 'd', 'B': '', 'C': 'f'}])
    >>> df2 = df1.copy()
    >>> df3 = df1.copy()
    >>> df1
       A  B  C
    0  a  b  c
    1  d     f

    >>> df1 = join_string_columns(df1, 'AB', ['A', 'B', 'C'])
    >>> df1
       A  B  C     AB
    0  a  b  c  a|b|c
    1  d     f    d|f
    
    Delete join_cols after operation.
    >>> df2 = join_string_columns(df2, 'AB', ['A', 'B', 'C'], delete=True, delimiter=',')
    >>> df2
          AB
    0  a,b,c
    1    d,f
    
    Reuse an existing column, but delete all other join_cols.
    >>> df3 = join_string_columns(df3, 'A', ['A', 'B', 'C'], delete=True, delimiter=',')
    >>> df3
           A
    0  a,b,c
    1    d,f
    """
    notnull = lambda s: pd.isnull(s) or str(s) == ""
    notnull = lambda x: not pd.isnull(x) and x != "nan"
    df[column] = df[join_cols].apply(
        lambda row: delimiter.join(filter(notnull, row)), axis=1
    )
    # df[column] = df[join_cols].apply(lambda row: delimiter.join(filter(None, row)), axis=1)

    if delete:
        # don't delete result column if it is in the list
        if column in join_cols:
            join_cols.remove(column)
        df.drop(columns=join_cols, inplace=True)

    return df


def backfill(df, column, columns, delete=False):
    df[column] = df[columns].bfill(axis=1).iloc[:, 0]
    if delete:
        if column in columns:
            columns.remove(column)
        df.drop(columns=columns, inplace=True)

    return df


# def join_columns(row, columns, delimiter='|'):
#     """Joins non-empty columns within a row"""
#     items = [row[col] for col in columns if len(row[col]) > 0]
#     return delimiter.join(items)


def assign_curie(identifier, curie):
    """Assigns a compact URI (CURIE) to an identifier"""
    if identifier == "":
        return identifier
    if ":" in identifier:
        # remove exiting prefix
        identifier = identifier.split(":")[1]

    return curie + ":" + identifier
