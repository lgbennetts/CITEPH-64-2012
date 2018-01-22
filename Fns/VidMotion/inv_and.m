% <!DOCTYPE html>
% <html lang='en'>
% <head>
% <meta charset='utf-8'>
% <meta content='GitLab Community Edition' name='description'>
% <title>
% dymat/inv_and.m at master - Dany Dumont / dymat | 
% GitLab
% </title>
% <link href="/assets/favicon-baaa14bade1248aa6165e9d34e7d83c0.ico" rel="shortcut icon" type="image/vnd.microsoft.icon" />
% <link href="/assets/application-2ac58be704fef592dcfa8a14124e2688.css" media="all" rel="stylesheet" />
% <link href="/assets/print-1df3ea9b8ff148a6745321899e0cb213.css" media="print" rel="stylesheet" />
% <script src="/assets/application-5f3c67da81199dd3da676a30746cf17f.js"></script>
% <meta content="authenticity_token" name="csrf-param" />
% <meta content="AjxTd2Sc5N2bP96n+2717LOdhj6JhDQkbzi0GkMrJOc=" name="csrf-token" />
% <script type="text/javascript">
% //<![CDATA[
% window.gon={};gon.default_issues_tracker="gitlab";gon.api_version="v3";gon.relative_url_root="";gon.default_avatar_url="https://gitlasso.uqar.ca/assets/no_avatar-adffbfe10d45b20495cd2a9b88974150.png";gon.current_user_id=3;gon.api_token="uQCzoNQqxppUxtS89RMy";
% //]]>
% </script>
% <meta content='width=device-width, initial-scale=1.0' name='viewport'>
% 
% 
% </head>
% 
% <body class='ui_modern dark_theme project' data-page='projects:blob:show' data-project-id='36'>
% 
% <header class='navbar navbar-fixed-top navbar-gitlab'>
% <div class='navbar-inner'>
% <div class='container'>
% <div class='app_logo'>
% <a class="home has_bottom_tooltip" href="/" title="Dashboard"><h1>GITLAB</h1>
% </a></div>
% <h1 class='title'><span><a href="/u/dumoda01">Dany Dumont</a> / <a href="/dumoda01/dymat">dymat</a></span></h1>
% <button class='navbar-toggle' data-target='.navbar-collapse' data-toggle='collapse' type='button'>
% <span class='sr-only'>Toggle navigation</span>
% <i class='fa fa-bars'></i>
% </button>
% <div class='navbar-collapse collapse'>
% <ul class='nav navbar-nav'>
% <li class='hidden-sm hidden-xs'>
% <div class='search'>
% <form accept-charset="UTF-8" action="/search" class="navbar-form pull-left" method="get"><div style="display:none"><input name="utf8" type="hidden" value="&#x2713;" /></div>
% <input class="search-input" id="search" name="search" placeholder="Search in this project" type="search" />
% <input id="group_id" name="group_id" type="hidden" />
% <input id="project_id" name="project_id" type="hidden" value="36" />
% <input id="search_code" name="search_code" type="hidden" value="true" />
% <input id="repository_ref" name="repository_ref" type="hidden" value="master" />
% 
% <div class='search-autocomplete-opts hide' data-autocomplete-path='/search/autocomplete' data-autocomplete-project-id='36' data-autocomplete-project-ref='master'></div>
% </form>
% 
% </div>
% <script>
%   $('.search-input').on('keyup', function(e) {
%     if (e.keyCode == 27) {
%       $('.search-input').blur()
%     }
%   })
% </script>
% 
% </li>
% <li class='visible-sm visible-xs'>
% <a class="has_bottom_tooltip" data-original-title="Search area" href="/search" title="Search"><i class='fa fa-search'></i>
% </a></li>
% <li>
% <a class="has_bottom_tooltip" data-original-title="Help" href="/help" title="Help"><i class='fa fa-question-circle'></i>
% </a></li>
% <li>
% <a class="has_bottom_tooltip" data-original-title="Public area" href="/explore" title="Explore"><i class='fa fa-globe'></i>
% </a></li>
% <li>
% <a class="has_bottom_tooltip" data-original-title="My snippets" href="/s/dumoda01" title="My snippets"><i class='fa fa-clipboard'></i>
% </a></li>
% <li>
% <a class="has_bottom_tooltip" data-original-title="New project" href="/projects/new" title="New project"><i class='fa fa-plus'></i>
% </a></li>
% <li>
% <a class="has_bottom_tooltip" data-original-title="Profile settings&quot;" href="/profile" title="Profile settings"><i class='fa fa-user'></i>
% </a></li>
% <li>
% <a class="has_bottom_tooltip" data-method="delete" data-original-title="Logout" href="/users/sign_out" rel="nofollow" title="Logout"><i class='fa fa-sign-out'></i>
% </a></li>
% <li class='hidden-xs'>
% <a class="profile-pic" href="/u/dumoda01" id="profile-pic"><img alt="User activity" src="https://gitlasso.uqar.ca//uploads/user/avatar/3/dany_cvismer_4.jpg" />
% </a></li>
% </ul>
% </div>
% </div>
% </div>
% </header>
% 
% 
% <script>
%   GitLab.GfmAutoComplete.dataSource = "/dumoda01/dymat/autocomplete_sources?type=NilClass&type_id=master%2Finv_and.m"
%   GitLab.GfmAutoComplete.setup();
% </script>
% 
% <div class='page-with-sidebar'>
% <div class='sidebar-wrapper'>
% <ul class='project-navigation nav nav-sidebar'>
% <li class="home"><a class="shortcuts-project" href="/dumoda01/dymat" title="Project"><i class='fa fa-dashboard'></i>
% <span>
% Project
% </span>
% </a></li><li class="active"><a class="shortcuts-tree" href="/dumoda01/dymat/tree/master"><i class='fa fa-files-o'></i>
% <span>
% Files
% </span>
% </a></li><li class=""><a class="shortcuts-commits" href="/dumoda01/dymat/commits/master"><i class='fa fa-history'></i>
% <span>
% Commits
% </span>
% </a></li><li class=""><a class="shortcuts-network" href="/dumoda01/dymat/network/master"><i class='fa fa-code-fork'></i>
% <span>
% Network
% </span>
% </a></li><li class=""><a class="shortcuts-graphs" href="/dumoda01/dymat/graphs/master"><i class='fa fa-area-chart'></i>
% <span>
% Graphs
% </span>
% </a></li><li class=""><a class="shortcuts-issues" href="/dumoda01/dymat/issues"><i class='fa fa-exclamation-circle'></i>
% <span>
% Issues
% <span class='count issue_counter'>0</span>
% </span>
% </a></li><li class=""><a class="shortcuts-merge_requests" href="/dumoda01/dymat/merge_requests"><i class='fa fa-tasks'></i>
% <span>
% Merge Requests
% <span class='count merge_counter'>0</span>
% </span>
% </a></li><li class=""><a class="shortcuts-wiki" href="/dumoda01/dymat/wikis/home"><i class='fa fa-book'></i>
% <span>
% Wiki
% </span>
% </a></li><li class="separate-item"><a class="stat-tab tab no-highlight" href="/dumoda01/dymat/edit"><i class='fa fa-cogs'></i>
% <span>
% Settings
% <i class='fa fa-angle-down'></i>
% </span>
% </a></li></ul>
% 
% </div>
% <div class='content-wrapper'>
% <div class='container-fluid'>
% <div class='content'>
% <div class='flash-container'>
% </div>
% 
% <div class='clearfix'>
% <div class='tree-ref-holder'>
% <form accept-charset="UTF-8" action="/dumoda01/dymat/refs/switch" class="project-refs-form" method="get"><div style="display:none"><input name="utf8" type="hidden" value="&#x2713;" /></div>
% <select class="project-refs-select select2 select2-sm" id="ref" name="ref"><optgroup label="Branches"><option selected="selected" value="master">master</option></optgroup><optgroup label="Tags"></optgroup></select>
% <input id="destination" name="destination" type="hidden" value="blob" />
% <input id="path" name="path" type="hidden" value="inv_and.m" />
% </form>
% 
% 
% </div>
% <div class='tree-holder' id='tree-holder'>
% <ul class='breadcrumb repo-breadcrumb'>
% <li>
% <i class='fa fa-angle-right'></i>
% <a href="/dumoda01/dymat/tree/master">dymat
% </a></li>
% <li>
% <a href="/dumoda01/dymat/blob/master/inv_and.m"><strong>
% inv_and.m
% </strong>
% </a></li>
% </ul>
% <ul class='blob-commit-info bs-callout bs-callout-info hidden-xs'>
% <li class='commit js-toggle-container'>
% <div class='commit-row-title'>
% <a class="commit_short_id" href="/dumoda01/dymat/commit/a47994590d328da7e8849ddcad7a54e1a8c7bdd2">a4799459</a>
% &nbsp;
% <span class='str-truncated'>
% <a class="commit-row-message" href="/dumoda01/dymat/commit/a47994590d328da7e8849ddcad7a54e1a8c7bdd2">Ajout de inv_and.m</a>
% </span>
% <a class="pull-right" href="/dumoda01/dymat/tree/a47994590d328da7e8849ddcad7a54e1a8c7bdd2">Browse Code »</a>
% <div class='notes_count'>
% </div>
% </div>
% <div class='commit-row-info'>
% <a class="commit-author-link has_tooltip" data-original-title="dany_dumont@uqar.ca" href="/u/dumoda01"><img alt="" class="avatar s16" src="https://gitlasso.uqar.ca//uploads/user/avatar/3/dany_cvismer_4.jpg" width="16" /> <span class="commit-author-name">Dany Dumont</span></a>
% <div class='committed_ago'>
% <time class='time_ago' data-placement='top' data-toggle='tooltip' datetime='2018-01-16T12:39:50Z' title='Jan 16, 2018 7:39am'>2018-01-16 07:39:50 -0500</time>
% <script>$('.time_ago').timeago().tooltip()</script>
%  &nbsp;
% </div>
% </div>
% </li>
% 
% </ul>
% <div class='tree-content-holder' id='tree-content-holder'>
% <article class='file-holder'>
% <div class='file-title clearfix'>
% <i class='fa fa-file'></i>
% <span class='file_name'>
% inv_and.m
% <small>334 Bytes</small>
% </span>
% <span class='options hidden-xs'><div class='btn-group tree-btn-group'>
% <a class="btn btn-small" href="/dumoda01/dymat/edit/master/inv_and.m">Edit</a>
% <a class="btn btn-small" href="/dumoda01/dymat/raw/master/inv_and.m" target="_blank">Raw</a>
% <a class="btn btn-small" href="/dumoda01/dymat/blame/master/inv_and.m">Blame</a>
% <a class="btn btn-small" href="/dumoda01/dymat/commits/master/inv_and.m">History</a>
% <a class="btn btn-small" href="/dumoda01/dymat/blob/a47994590d328da7e8849ddcad7a54e1a8c7bdd2/inv_and.m">Permalink</a>
% </div>
% <button class="remove-blob btn btn-small btn-remove" data-target="#modal-remove-blob" data-toggle="modal" name="button" type="submit">Remove
% </button></span>
% </div>
% <div class='file-content code'>
% <div class='dark highlighted-data'>
% <div class='line-numbers'>
% <a href="#L1" id="L1" rel="#L1"><i class='fa fa-link'></i>
% 1
% </a><a href="#L2" id="L2" rel="#L2"><i class='fa fa-link'></i>
% 2
% </a><a href="#L3" id="L3" rel="#L3"><i class='fa fa-link'></i>
% 3
% </a><a href="#L4" id="L4" rel="#L4"><i class='fa fa-link'></i>
% 4
% </a><a href="#L5" id="L5" rel="#L5"><i class='fa fa-link'></i>
% 5
% </a><a href="#L6" id="L6" rel="#L6"><i class='fa fa-link'></i>
% 6
% </a><a href="#L7" id="L7" rel="#L7"><i class='fa fa-link'></i>
% 7
% </a><a href="#L8" id="L8" rel="#L8"><i class='fa fa-link'></i>
% 8
% </a><a href="#L9" id="L9" rel="#L9"><i class='fa fa-link'></i>
% 9
% </a><a href="#L10" id="L10" rel="#L10"><i class='fa fa-link'></i>
% 10
% </a><a href="#L11" id="L11" rel="#L11"><i class='fa fa-link'></i>
% 11
% </a><a href="#L12" id="L12" rel="#L12"><i class='fa fa-link'></i>
% 12
% </a><a href="#L13" id="L13" rel="#L13"><i class='fa fa-link'></i>
% 13
% </a><a href="#L14" id="L14" rel="#L14"><i class='fa fa-link'></i>
% 14
% </a><a href="#L15" id="L15" rel="#L15"><i class='fa fa-link'></i>
% 15
% </a><a href="#L16" id="L16" rel="#L16"><i class='fa fa-link'></i>
% 16
% </a><a href="#L17" id="L17" rel="#L17"><i class='fa fa-link'></i>
% 17
% </a><a href="#L18" id="L18" rel="#L18"><i class='fa fa-link'></i>
% 18
% </a><a href="#L19" id="L19" rel="#L19"><i class='fa fa-link'></i>
% 19
% </a><a href="#L20" id="L20" rel="#L20"><i class='fa fa-link'></i>
% 20
% </a><a href="#L21" id="L21" rel="#L21"><i class='fa fa-link'></i>
% 21
% </a></div>
% <div class='highlight'>
% <pre><code class='inv_and.m'>

function v = inv_and(u)
% INV_AND - This function solves the linear problem 
% u must be a 1x4 or 4x1 vector

if size(u) == [4 1]
elseif size(u) == [1 4];
    u = u';
else
    disp([' u must be a 1x4 or 4x1 vector'])
    return
end

M = [1 0 0 0; ...
     1 1 0 0; ...
     0 1 1 0; ...
     0 0 1 1; ...
     0 0 0 1];

v = M*u;

end

% </code></pre>
% </div>
% </div>
% 
% </div>
% 
% </article>
% </div>
% 
% </div>
% <div class='modal hide' id='modal-remove-blob'>
% <div class='modal-dialog'>
% <div class='modal-content'>
% <div class='modal-header'>
% <a class='close' data-dismiss='modal' href='#'>×</a>
% <h3 class='page-title'>Remove inv_and.m</h3>
% <p class='light'>
% From branch
% <strong>master</strong>
% </p>
% </div>
% <div class='modal-body'>
% <form accept-charset="UTF-8" action="/dumoda01/dymat/blob/master/inv_and.m" class="form-horizontal" method="post"><div style="display:none"><input name="utf8" type="hidden" value="&#x2713;" /><input name="_method" type="hidden" value="delete" /><input name="authenticity_token" type="hidden" value="AjxTd2Sc5N2bP96n+2717LOdhj6JhDQkbzi0GkMrJOc=" /></div>
% <div class='form-group commit_message-group'>
% <label class="control-label" for="commit_message">Commit message
% </label><div class='col-sm-10'>
% <div class='commit-message-container'>
% <div class='max-width-marker'></div>
% <textarea class="form-control" id="commit_message" name="commit_message" placeholder="Removed this file because..." required="required" rows="3">
% </textarea>
% </div>
% </div>
% </div>
% 
% <div class='form-group'>
% <div class='col-sm-2'></div>
% <div class='col-sm-10'>
% <button class="btn btn-remove btn-remove-file" name="button" type="submit">Remove file</button>
% <a class="btn btn-cancel" data-dismiss="modal" href="#">Cancel</a>
% </div>
% </div>
% </form>
% 
% </div>
% </div>
% </div>
% </div>
% <script>
%   disableButtonIfEmptyField('#commit_message', '.btn-remove-file')
% </script>
% 
% 
% </div>
% </div>
% </div>
% </div>
% </div>
% 
% 
% </body>
% </html>
