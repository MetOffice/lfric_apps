# PR Summary

Sci/Tech Reviewer: <!-- SR id, filled by author when ready for review (e.g. @octocat) -->
Code Reviewer: <!-- CR id, filled by SSD/CCD (e.g. @octocat) -->

<!-- To be completed by the developer -->

<!-- Provide a brief description of the changes in this PR, including any notes
     useful for reviewers -->

<!-- List any linked PRs here
- linked MetOffice/<REPO-NAME>#<pr-number>
-->

<!-- List any blocking PRs or issues to be closed here
- is blocked-by #pr-number
- blocks #pr-number
- closes #issue-number (auto-closes the issue)
- fixes #issue-number (auto-closes the issue)
- is related to #issue-number
-->

## Code Quality Checklist

- [ ] I have performed a self-review of my own code
- [ ] My code follows the project's [style guidelines](https://metoffice.github.io/lfric_core/how_to_contribute/index.html#how-to-contribute-index)
- [ ] Comments have been included that aid understanding and enhance the readability of the code
- [ ] My changes generate no new warnings
- [ ] All automated checks in the CI pipeline have completed successfully

## Testing

- [ ] I have tested this change locally, using the LFRic Apps rose-stem suite
- [ ] If any tests fail (rose-stem or CI) the reason is understood and acceptable (e.g. kgo changes)
- [ ] I have added tests to cover new functionality as appropriate (e.g. system tests, unit tests, etc.)
- [ ] Any new tests have been assigned an appropriate amount of compute resource and have been allocated to an appropriate testing group (i.e. the developer tests are for jobs which use a small amount of compute resource and complete in a matter of minutes)

<!-- Describe other testing performed (if applicable) -->

### trac.log

<!-- Paste your trac.log from testing output here -->

## Security Considerations

- [ ] I have reviewed my changes for potential security issues
- [ ] Sensitive data is properly handled (if applicable)
- [ ] Authentication and authorisation are properly implemented (if applicable)

## Performance Impact

- [ ] Performance of the code has been considered and, if applicable, suitable performance measurements have been conducted

## AI Assistance and Attribution

- [ ] Some of the content of this change has been produced with the assistance of _Generative AI tool name_ (e.g., Met Office Github Copilot Enterprise, Github Copilot Personal, ChatGPT GPT-4, etc) and I have followed the [Simulation Systems AI policy](https://metoffice.github.io/simulation-systems/FurtherDetails/ai.html) (including attribution labels)

<!-- If AI has been used, please provide more details here -->

## Documentation

- [ ] Where appropriate I have updated documentation related to this change and confirmed that it builds correctly

## PSyclone Approval

- [ ] If you have edited any PSyclone-related code (e.g. PSyKAl-lite, Kernel interface, optimisation scripts, LFRic data structure code) then please contact the [TCD Team](mailto:ToolsCollabDevTeam@metoffice.gov.uk)

# Sci/Tech Review

<!-- To be completed by the Sci/Tech Reviewer -->
<!-- May be skipped for trivial tickets -->

- [ ] I understand this area of code and the changes being added
- [ ] The proposed changes correspond to the pull request description
- [ ] Documentation is sufficient (do documentation papers need updating)
- [ ] Sufficient testing has been completed

(_Please alert the code reviewer via a tag when you have approved the SR_)

# Code Review

<!-- To be completed by the Code Reviewer -->

- [ ] All dependencies have been resolved
- [ ] Related Issues have been properly linked and addressed
- [ ] CLA compliance has been confirmed
- [ ] Code quality standards have been met
- [ ] Tests are adequate and have passed
- [ ] Documentation is complete and accurate
- [ ] Security considerations have been addressed
- [ ] Performance impact is acceptable
