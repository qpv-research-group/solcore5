import inspect

bibliography: list = []
bibliography_id: list = []
science_tracking: list = []
track_science_references: bool = False
track_each_reference_call: bool = False


def hexID(obj: str) -> str:
    return "{}".format(id(obj))


def science_reference(short_purpose: str, reference: str) -> None:
    """Marker acting as a reference for the origin of specific information

    Acts as a marker in code where particular alogrithms/data/... originates. General execution of code silently passes
    these markers, but remembers how and where they were called. Which markers were passed in a particular program run
    can be recalled with print_references().

    Arguments:
    short_purpose: Identify the thing being referenced (string)
    reference: The reference itself, in any sensible format.
    """
    global bibliography, bibliography_id, science_tracking, track_science_references
    if not track_science_references:
        return
    stack = inspect.stack()
    frame, path, line, function, context, index = stack[1]
    del stack

    arguments = inspect.formatargvalues(*inspect.getargvalues(frame))

    availableObjects = {}
    availableObjects.update(frame.f_globals)
    availableObjects.update(frame.f_locals)
    del frame

    func_id = hexID(availableObjects[function])
    identifier = "{} [{}]".format(function, func_id)
    call_record = "{}{}".format(function, arguments)

    add_anyway = False
    if identifier not in bibliography_id:
        bibliography_id.append(identifier)
        bibliography.append(reference)
        add_anyway = True
    if track_each_reference_call or add_anyway:
        science_tracking.append((call_record, short_purpose, bibliography_id.index(identifier) + 1))


def print_references() -> None:
    """ recall the science_reference markers passed, print out the referenes."""
    global bibliography, bibliography_id, science_tracking, track_each_reference_call

    "List of references encountered while executing"
    if not track_each_reference_call: print(
        "showing first call only - override with solcore.track_science(track_each_call=True)")

    for record, purpose, index in science_tracking:
        print("[{}] {} - {}".format(index, purpose, record))
    print()
    for i, b in enumerate(bibliography):
        print("[{}] - {}".format(i + 1, b))


def track_science(track_each_call: bool = False):
    """configure science references -- determine whether or not each call separately or only the first for each
    reference."""
    global track_science_references, track_each_reference_call
    track_science_references = True
    track_each_reference_call = track_each_call


if __name__ == "__main__":
    def myRandomFunctionA() -> None:
        myRandomFunctionB(d=1)


    def myRandomFunctionB(d) -> None:
        a = 1

        science_reference("doing something smart", "My Awesome Book, by me.")


    track_science()
    myRandomFunctionA()
    myRandomFunctionA()
    print_references()
