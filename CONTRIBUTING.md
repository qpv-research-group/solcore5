# Contributing to Solcore

Thanks for taking the time to contribute to Solcore!

The following is a set of guidelines for contributing to Solcore, a multi-scale, Python-based library for the modelling of solar cells and semiconductor materials. The goal of these guidelines is to make the development of the package efficient and sustainable and to ensure that every commit makes it better, more readable, more robust and better documented. Please, feel free suggest changes and improvements. 

(this guide is based on the [Atom editor guide](https://github.com/atom/atom/blob/master/CONTRIBUTING.md))

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Your First Code Contribution](#your-first-code-contribution)
  * [Pull Requests](#pull-requests)

[Styleguides](#styleguides)
  * [Git Commit Messages](#git-commit-messages)
  * [Documentation Styleguide](#documentation-styleguide)

## Code of Conduct

This project and everyone participating in it is governed by the [Solcore Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [the repository Administrator](mailto:nekins@unsw.edu.au).

## How Can I Contribute?

### Reporting Bugs

This section guides you through submitting a bug report for Solcore. Following these guidelines helps maintainers and the community understand your report :pencil:, reproduce the behavior :computer: :computer:, and find related reports :mag_right:.

Before creating bug reports, please check [this list](https://github.com/dalonsoa/solcore5/issues) (including the closed issues) as you might find out that you don't need to create one. When you are creating a bug report, please [include as many details as possible](#how-do-i-submit-a-good-bug-report). 

> **Note:** If you find a **Closed** issue that seems like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one.

#### How Do I Submit A (Good) Bug Report?

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/). Explain the problem and include additional details to help maintainers reproduce the problem:

* **Use a clear and descriptive title** for the issue to identify the problem.
* **Describe the exact steps which reproduce the problem** in as many details as possible. For example, start by explaining how you installed Solcore and what you where trying to do.
* **Provide specific examples to demonstrate the steps**. Include links to files or GitHub projects, or copy/pasteable snippets, which you use in those examples. If you're providing snippets in the issue, use [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the behavior you observed after following the steps** and point out what exactly is the problem with that behavior.
* **Explain which behavior you expected to see instead and why.**
* **If there is any error output in the temrinal, include that output with your report.**

Provide more context by answering these questions:

* **Did the problem start happening recently** (e.g. after updating to a new version of Solcore) or was this always a problem?
* If the problem started happening recently, **can you reproduce the problem in an older version of Solcore?** What's the most recent version in which the problem doesn't happen? You can download older versions of Solcore from [the releases page](https://github.com/dalonsoa/solcore5/releases).
* **Can you reliably reproduce the issue?** If not, provide details about how often the problem happens and under which conditions it normally happens.

Include details about your configuration and environment:

* **Which version of Solcore are you using?** 
* **What's the name and version of the OS you're using**?
* **Are you running Solcore in a virtual machine?** If so, which VM software are you using and which operating systems and versions are used for the host and the guest?
* **If the problem is related with the Poisson-drift-diffusion solver**, which compiler did you use? What were the installation steps you followed?

### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for Solcore, including completely new features and minor improvements to existing functionality. Following these guidelines helps maintainers and the community understand your suggestion :pencil: and find related suggestions :mag_right:.

Before creating enhancement suggestions, please check [this list](https://github.com/dalonsoa/solcore5/issues) (including closed issues) as you might find out that you don't need to create one. When you are creating an enhancement suggestion, please [include as many details as possible](#how-do-i-submit-a-good-enhancement-suggestion).

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/). Create an issue on that repository and provide the following information:

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Explain why this enhancement would be useful** to most Solcore users, maybe including some links to scientific papers showing the enhancement in action.
* **List some other packages or applications where this enhancement exists.**
* **Specify the name and version of the OS you're using.**

### Your First Code Contribution

Unsure where to begin contributing to Atom? You can start by looking through these `beginner` and `help-wanted` issues:

* [Beginner issues](https://github.com/dalonsoa/solcore5/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22) - issues which should only require a few lines of code, and a test or two.
* [Help wanted issues](https://github.com/dalonsoa/solcore5/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22) - issues which should be a bit more involved than `beginner` issues.

### Pull Requests

The process described here has several goals:

- Maintain Solcore's quality
- Fix problems that are important to users
- Engage the community in working toward the best possible Solcore
- Enable a sustainable system for Solcore's maintainers to review contributions

Please follow these steps to have your contribution considered by the maintainers:

1. **Describe clearly what is the purpose of the pull request**. Refer to the relevant issues on [Bugs](#reporting-bugs) or [Enhancements](#suggesting-enhancements). In general, an issue should always be open *prior* to a pull request, to discuss its contents with a maintainer and make sure it makes sense for Solcore. If the pull request is a work in progress that will take some time to be ready but still you want to discuss it with the community, open a [draft pull request](https://github.blog/2019-02-14-introducing-draft-pull-requests/). 
2. **Include relevant unit tests and integration tests, where needed**. Solcore's test suite is quite limited at the moment. We are working to improve this and tests as many features as possible, so any new addition to the code must come with its own set of tests to avoid going backwards in this matter.   
3. **For new features and enhancements, include documentation and examples**. Both in the code, as docstrings in classes, functions and modules, and as proper documentation describing how to use the new feature. 
4. Follow the [styleguides](#styleguides)
5. After you submit your pull request, verify that all [status checks](https://help.github.com/articles/about-status-checks/) are passing <details><summary>What if the status checks are failing?</summary>If a status check is failing, and you believe that the failure is unrelated to your change, please leave a comment on the pull request explaining why you believe the failure is unrelated. A maintainer will re-run the status check for you. If we conclude that the failure was a false positive, then we will open an issue to track that problem with our status check suite.</details>

While the prerequisites above must be satisfied prior to having your pull request reviewed, the reviewer(s) may ask you to complete additional design work, tests, or other changes before your pull request can be ultimately accepted.

## Styleguides

### Git Commit Messages

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include `[ci skip]` in the commit title
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :non-potable_water: `:non-potable_water:` when plugging memory leaks
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :lock: `:lock:` when dealing with security
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies
    * :shirt: `:shirt:` when removing linter warnings
    
### Documentation Styleguide

* Use [Markdown](https://daringfireball.net/projects/markdown).
* Reference methods and classes in markdown with the custom `{}` notation:
    * Reference classes with `{ClassName}`
    * Reference instance methods with `{ClassName::methodName}`
    * Reference class methods with `{ClassName.methodName}`
