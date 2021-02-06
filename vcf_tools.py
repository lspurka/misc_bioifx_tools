# tools for handling VCFs


import os


def check_and_reformat_user_inputs(vcf_path, vcf_fields, start_position, end_position, pad, ref_alleles, alt_alleles, 
        qual_threshold, variant_types):
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
        vcf_fields = [i.upper() for i in vcf_fields]
    
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
    
    # TODO remove
    print("YAY #1")

    return vcf_fields, ref_alleles, alt_alleles, variant_types


def filter_vcf_df(vcf_df, vcf_fields, start_position, end_position, pad, ref_alleles, alt_alleles, qual_threshold, 
        variant_types):
    """ Filter VCF dataframe """

    if start_position:
        assert type(start_position)
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
        for col in vcf_df.columns:
            if "ALT" in col:
                vcf_df = vcf_df.query(f"{col} in @alt_alleles")

    if qual_threshold:
        vcf_df = vcf_df.query("QUAL >= @qual_threshold")

    if variant_types:
        df = pd.DataFrame()
        if "SUBSTITUTION" in variant_types:
            filtered_df = vcf_df.query("POS == END")
            df = pd.concat([df, filtered_df])
        if ("DELETION" or "DELINS") in variant_types:


def add_var_seq_is_diff_col()
    # TODO left off


def get_vcf_df(vcf_path, vcf_fields):
    """ Returns a non-empty dataframe of the parsed VCF 
    :param str vcf_path: path to VCF file
    :param list vcf_fields: None or list of VCF fields to output
    :return pandas.DataFrame: dataframe of VCF variants
    """

    # Check that VCF file has variants
    # if vcf_fields is set to None, the default columns from allel will be output
    vcf_df = allel.vcf_to_dataframe(vcf_path, fields=vcf_fields)
    assert not vcf_df.empty, f"ERROR: No VCF contents for: {vcf_path}"

    # Check that columns from parsed VCF have variant rows
    vcf_df.dropna(axis=1, how="all", inplace=True)
    if vcf_fields:
        assert not vcf_df.empty, "ERROR: No results, Double check the specified vcf_fields."
    assert not vcf_df.empty, "ERROR: No results.  Check inputs."

    vcf_df.columns = vcf_df.columns.str.upper()

    # Need REF column for downstream filtering
    assert "REF" in vcf_df.columns, "ERROR: Need REF column."

    return vcf_df


def get_vcf_variants(vcf_path, vcf_fields=None, start_position=None, end_position=None, pad=None, ref_alleles=None, 
        alt_alleles=None, qual_threshold=None, variant_types=None):
    """ Tool to get variants in a VCF based on any of: start position, end position, ref allele, alt allele, quality 
    score. If no filters are specified, the function will return a dataframe of the parsed VCF.  This tool works best 
    for small scale variants, rather than CNV VCFs.
    :param str vcf_path: path to VCF file
    :param list vcf_fields: optional, list of VCF fields to output
    :param int start_position: optional, start position to filter variants by.
    :param int end_position: optional, end position to filter variants by.
    :param int pad: optional, pads the pos and end param if either are specified.
    :param list(str) ref_alleles: optional, specific reference alleles to filter variants by
    :param list(str) alt_alleles: optional, specific alternate alleles to filter variants by
    :param float qual_threshold: optional, quality score threshold, return results at or above the threshold.
    :param list(str) variant_types: optional, list of variant types to filter variants by, any of: 'substitution', 
        'indel', 'deletion', 'insertion', 'delins'. Delins variants will be included if any of 'indel', 'deletion', 
        'insertion', 'delins' are specified.
    :return pandas.DataFrame: dataframe of VCF variants
    """


    # Check user input, reformat lists to uppercase
    vcf_fields, ref_alleles, alt_alleles, variant_types = _check_and_reformat_user_inputs(vcf_path=vcf_path, 
            vcf_fields=vcf_fields, start_position=start_position, end_position=end_position, pad=pad, 
            ref_alleles=ref_alleles, alt_alleles=alt_alleles, qual_threshold=qual_threshold, 
            variant_types=variant_types)

    # Get non-empty dataframe from parsed VCF
    vcf_df = get_vcf_df(vcf_path=vcf_path, vcf_fields=vcf_fields)

    # Add END position column
    vcf_df["END"] = vcf_df["POS"] + vcf_df["REF"].str.len() - 1 # subtract one since VCFs are 1-indexed

    # Add sequence length columns REF_LEN and ALT_LEN, which are the lengths of the allele sequences
    vcf_df["REF_LEN"] = vcf_df["REF"].str.len()
    vcf_df["ALT_LEN"] = vcf_df["ALT"].str.len()

    # Add variant sequence is different from reference column (will be helpful for determining delins)
    vcf_df = add_var_seq_is_diff_col(vcf_df)


    filtered_df = _filter_vcf_df(vcf_df=vcf_df, start_position=start_position, end_position=end_position, pad=pad, 
            ref_alleles=ref_alleles, alt_alleles=alt_alleles, qual_threshold=qual_threshold, 
            variant_types=variant_types)








# TODO tests


