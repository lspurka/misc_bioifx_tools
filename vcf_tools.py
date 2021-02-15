# tools for handling VCFs

import allel
import os


def get_vcf_df(vcf_path, vcf_fields=None, alt_number=1):
    """ Returns a non-empty dataframe of the parsed VCF 
    :param str vcf_path: path to VCF file
    :param list vcf_fields: optional, list of additional VCF fields to output, default is:
        ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER_PASS']
    :param int alt_number: Looks for the number of specified alt alleles in the VCF
    :return pandas.DataFrame: dataframe of VCF variants
    """

    # Check that VCF file has variants
    # if vcf_fields is set to None, the default columns from allel will be output

    # TODO fix below logic to be more understandable
    default_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER_PASS"]
    if vcf_fields == ["*"]:
        vcf_fields = "*"
    elif vcf_fields:
        vcf_fields = default_fields + vcf_fields
    else:
        vcf_fields = default_fields

    vcf_df = allel.vcf_to_dataframe(vcf_path, fields=vcf_fields, alt_number=alt_number)
    assert not vcf_df.empty, f"ERROR: No VCF contents for: {vcf_path}"

    # Check that columns from parsed VCF have variant rows
    # vcf_df.dropna(axis=1, how="all", inplace=True)
    # if vcf_fields:
    #     assert not vcf_df.empty, "ERROR: No results, Double check the specified vcf_fields."
    # assert not vcf_df.empty, "ERROR: No results.  Check inputs."

    return vcf_df


def add_columns_to_vcf_df(vcf_df):
    """ TODO docstring """

    def _get_var_type(row):
        """ TODO add docstring """
        # TODO double-check and test below:
        if row["REF_LEN"] == row["ALT_LEN"] == 1:
            return "snv"
        elif row["REF_LEN"] > 1 and row["ALT_LEN"] > 1:
            return "delins"
        elif row["REF_LEN"] > 1:
            return "deletion"
        elif row["ALT_LEN"] > 1:
            return "insertion"
        return "N/A"

    # Add allele sequence length columns for the ref and alt allele(s)
    vcf_df["REF_LEN"] = vcf_df["REF"].str.len()
    vcf_df["ALT_LEN"] = vcf_df["ALT"].str.len()

    # Add END position column
    vcf_df["END"] = vcf_df["POS"] + vcf_df["REF_LEN"] - 1 # subtract one since VCFs are 1-indexed

    # Add variant type column, which determines if the variant is specifically one of: 'snv', 'delins', 'deletion', 
    # 'insertion'
    vcf_df["VAR_TYPE"] = vcf_df.apply(_get_var_type, axis=1)

    return vcf_df


def filter_vcf_df(vcf_df, chrom, start_position, end_position, pad, ref_alleles, alt_alleles, qual_threshold, 
        variant_types):
    """ Filter VCF dataframe 
    # TODO finish docstring
    """

    if chrom:
        vcf_df = vcf_df.query("CHROM == @chrom")

    if start_position:
        if pad:
            start_position -= pad
        vcf_df = vcf_df.query("POS >= @start_position")
    
    if end_position:
        if pad:
            end_position += pad
        vcf_df = vcf_df.query("END <= @end_position")

    if ref_alleles:
        vcf_df = vcf_df.query("REF in @ref_alleles")

    if alt_alleles:
        vcf_df = vcf_df.query("ALT in @alt_alleles")

    if qual_threshold:
        vcf_df = vcf_df.query("QUAL >= @qual_threshold")

    if variant_types:
        var_types_to_use = []
        # TODO double-check that this is correct:
        if "snv" in variant_types:
            var_types_to_use.append("snv")
        if "deletion" or "indel" or "delins" in variant_types:
            var_types_to_use.append("deletion")
        if "insertion" or "indel" or "delins" in variant_types:
            var_types_to_use.append("insertion")
        if "delins" in variant_types:
            var_types_to_use.append("delins")
        var_types_to_use = list(set(var_types_to_use))
        vcf_df.query("VAR_TYPE in @var_types_to_use")

    return vcf_df


def get_vcf_variants(vcf_path, vcf_fields=None, alt_number=1, chrom=None, start_position=None, end_position=None, 
        pad=None, ref_alleles=None, alt_alleles=None, qual_threshold=None, variant_types=None):
    """ Tool to get variants in a VCF based on any of: start position, end position, ref allele, alt allele, quality 
    score, variant type. If no filters are specified, the function will return a dataframe of the parsed VCF.  This 
    tool works best for small scale variants, rather than CNV VCFs.
    :param str vcf_path: path to VCF file
    :param list vcf_fields: optional, list of additional VCF fields to output, default is:
        ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER_PASS'].  If the user wants all columns output, use ["*"]
    # TODO add in 'alt_number' functionality. `numalt` is a helpful column output by allel.vcf_to_dataframe
    :param list alt_number: optional, if >1, ALT columns will be named 'ALT_1', 'ALT_2', 'ALT_3', etc, otherwise just 
        one will be named 'ALT'
    :param str chrom: optional, chromosome to filter variants by.
    :param int start_position: optional, start position to filter variants by.
    :param int end_position: optional, end position to filter variants by.
    :param int pad: optional, pads the pos and end param if either are specified.
    :param list(str) ref_alleles: optional, specific reference alleles to filter variants by
    :param list(str) alt_alleles: optional, specific alternate alleles to filter variants by
    :param float qual_threshold: optional, quality score threshold, return results at or above the threshold.
    :param list(str) variant_types: optional, list of variant types to filter variants by, any of: 'snv', 
        'indel', 'deletion', 'insertion', 'delins'. Delins variants will be included if any of 'indel', 'deletion', 
        'insertion', 'delins' are specified.
    :return pandas.DataFrame: dataframe of VCF variants
    """

    # TODO add custom_filters param, which is a dictionary where the key is the VCF column name (for example, "AB" or 
    # "interpretation") and has any of the following sub keys "eq" (equal to), "ne" (not equal to) "gt" (greater than), 
    # "lt" (less than), "gte" (greater than or equal to), "lte" (less than or equal to)
    # So for AB and interpretation, the input might look like: 
    # {"AB": {"gte":  0.1}, "interpretation": {"eq": ["pathogenic", "likely pathogenic"]}}
    # check that the input types, for example import numbers, isinstance(5.4534, numbers.Real).  eq / ne can be lists, 
    # numbers, or strings

    # Check user input, reformat lists to uppercase
    vcf_fields, chrom, alt_number, ref_alleles, alt_alleles, variant_types = _check_and_reformat_user_inputs(
            vcf_path=vcf_path, vcf_fields=vcf_fields, alt_number=alt_number, chrom=chrom, start_position=start_position, 
            end_position=end_position, pad=pad, ref_alleles=ref_alleles, alt_alleles=alt_alleles, 
            qual_threshold=qual_threshold, variant_types=variant_types)

    # Get non-empty dataframe from parsed VCF
    vcf_df = get_vcf_df(vcf_path=vcf_path, vcf_fields=vcf_fields, alt_number=alt_number)

    # Add columns to enable filtering
    vcf_df = add_columns_to_vcf_df(vcf_df)


    vcf_df = filter_vcf_df(vcf_df=vcf_df, chrom=chrom, start_position=start_position, end_position=end_position, 
            pad=pad, ref_alleles=ref_alleles, alt_alleles=alt_alleles, qual_threshold=qual_threshold, 
            variant_types=variant_types)
    
    return vcf_df


def _check_and_reformat_user_inputs(vcf_path, vcf_fields, chrom, alt_number, start_position, end_position, pad, 
        ref_alleles, alt_alleles, qual_threshold, variant_types):
    """ Check inputs to function get_vcf_variants """

    def _generate_error_message(param_str, type_str):
        """ Helper function to generate error message.
        :param str param_str: param name
        :param str type_str: type that param should be
        :return str: error message
        """
        return f"ERROR: '{param_str}' param needs to be {type_str}."


    # VCF checks
    # Check that VCF path exists
    assert os.path.exists(vcf_path), f"ERROR: 'vcf_path' does not exist: {vcf_path}"
    assert vcf_path.endswith(".vcf"), f"ERROR: 'vcf_path' does not have correct extension '.vcf': {vcf_path}"

    if vcf_fields:
        error_message = _generate_error_message("vcf_fields", "a list of strings")
        assert type(vcf_fields) == list, error_message
        assert set([type(i) for i in vcf_fields]) == {str}, error_message

    if chrom:
        chrom = str(chrom)

    if alt_number:
        assert type(alt_number) == int, _generate_error_message("alt_number", "an integer")
    
    if start_position:
        assert type(start_position) == int, _generate_error_message("start_position", "an integer")

    if end_position:
        assert type(end_position) == int, _generate_error_message("end_position", "an integer")
    
    if start_position and end_position:
        assert start_position <= end_position, "ERROR: start_position should be less than or equal to the " \
                "end_position."

    if pad:
        assert type(pad) == int, _generate_error_message("pad", "an integer")
    
    if ref_alleles:
        error_message = _generate_error_message("ref_alleles", "a list of strings")
        assert type(ref_alleles) == list, error_message
        assert set([type(i) for i in ref_alleles]) == {str}, error_message
        ref_alleles = [i.upper() for i in ref_alleles]

    if alt_alleles:
        error_message = _generate_error_message("alt_alleles", "a list of strings")
        assert type(alt_alleles) == list, error_message
        assert set([type(i) for i in alt_alleles]) == {str}, error_message
        alt_alleles = [i.upper() for i in alt_alleles]

    if qual_threshold:
        assert type(qual_threshold) in [float, int], _generate_error_message("qual_threshold", "a float or integer")

    if variant_types:
        error_message = _generate_error_message("variant_types", "a list of strings")
        assert type(variant_types) == list, error_message
        assert set([type(i) for i in variant_types]) == {str}, error_message
        variant_types = [i.upper() for i in variant_types]

    return vcf_fields, chrom, alt_number, ref_alleles, alt_alleles, variant_types



    


# TODO tests
# TODO add identifier for filtering in ID, for example for rsIDs
# TODO add outfile optional param
# TODO accommodate alt_number > 1



