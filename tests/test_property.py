def test_property_class():
    from solcore.parameter import Parameter

    answer = Parameter(42, description="The answer to everything")
    assert answer.magnitude == 42
    assert str(answer.units) == "dimensionless"
    assert answer.description == "The answer to everything"

    answer = Parameter(
        8848, "m", description="Height of mount Everest", reference="Wikipedia"
    )
    assert answer.magnitude == 8848
    assert str(answer.units) == "meter"
    assert answer.description == "Height of mount Everest"
    assert answer.reference == "Wikipedia"

    answer = Parameter("8848 m", description="Height of mount Everest")
    assert answer.magnitude == 8848
    assert str(answer.units) == "meter"
    assert answer.description == "Height of mount Everest"

    answer = Parameter("1 ns")
    assert answer.magnitude == 1
    assert str(answer.units) == "nanosecond"
    assert answer.description == ""
